//! Time-resolved magnetic/plasma providers and bounded local integration.
//! Fields and paths must share one Cartesian frame and coordinate time.
//! No geomagnetic, ionospheric or weather dataset is bundled or fetched.
use crate::frames::{Epoch, FrameGraph, FrameId, Rotation, Vec3};
use crate::{Error, Result, finite};

pub trait MagneticField: Send + Sync {
    fn teslas_at(&self, position_m: Vec3, epoch: Epoch) -> Result<Vec3>;
}
pub trait ElectronDensity: Send + Sync {
    fn per_m3_at(&self, position_m: Vec3, epoch: Epoch) -> Result<f64>;
}
/// Rotate a BODY/local vector field into graph root axes, sampling its local position.
/// A magnetic vector transforms by rotation, never by translation or frame velocity.
pub struct FramedMagneticField<'a> {
    pub graph: &'a FrameGraph,
    pub frame: FrameId,
    pub field: &'a dyn MagneticField,
}
impl MagneticField for FramedMagneticField<'_> {
    fn teslas_at(&self, position_m: Vec3, epoch: Epoch) -> Result<Vec3> {
        position_m.validate()?;
        let pose = self
            .graph
            .relative_pose(self.frame, self.graph.root(), epoch)?;
        let local = pose
            .rotation
            .inverse()
            .apply(position_m - pose.origin.position_m);
        let field = self.field.teslas_at(local, epoch)?;
        field.validate()?;
        let result = pose.rotation.apply(field);
        result.validate()?;
        Ok(result)
    }
}

/// Ideal static dipole; educational/local test model, not IGRF/WMM or a storm model.
/// Uses the conventional approximate SI coefficient mu0/(4 pi) = 1e-7.
pub struct DipoleField {
    pub centre_m: Vec3,
    pub moment_a_m2: Vec3,
    pub minimum_radius_m: f64,
}
impl MagneticField for DipoleField {
    fn teslas_at(&self, position_m: Vec3, _: Epoch) -> Result<Vec3> {
        position_m.validate()?;
        self.centre_m.validate()?;
        self.moment_a_m2.validate()?;
        finite(self.minimum_radius_m, "dipole coverage radius")?;
        if self.minimum_radius_m <= 0.0 {
            return Err(Error::InvalidInput(
                "dipole coverage radius must be positive".into(),
            ));
        }
        let delta = position_m - self.centre_m;
        let radius = delta.norm();
        finite(radius, "dipole radius")?;
        if radius < self.minimum_radius_m {
            return Err(Error::OutsideCoverage);
        }
        let n = delta * (1.0 / radius);
        let field = (n * (3.0 * self.moment_a_m2.dot(n)) - self.moment_a_m2)
            * (1e-7 / radius / radius / radius);
        field.validate()?;
        Ok(field)
    }
}

/// Local right-handed basis: z along B, x along the perpendicular projection of a
/// supplied guide. Unlike geographic longitude, this works at both geographic poles.
/// A zero B or parallel guide has no unique basis and returns an error.
pub fn field_aligned_basis(field_t: Vec3, guide: Vec3) -> Result<Rotation> {
    let z = field_t.normalized()?;
    let guide = guide.normalized()?;
    let perpendicular = guide - z * guide.dot(z);
    if perpendicular.norm() < 1e-12 {
        return Err(Error::InvalidInput(
            "field-aligned guide is parallel to B".into(),
        ));
    }
    let x = perpendicular.normalized()?;
    let y = z.cross(x);
    Rotation::new([[x.x, y.x, z.x], [x.y, y.y, z.y], [x.z, y.z, z.z]])
}

/// Frozen-time field line, integrated with RK4 in signed arc-length steps.
/// All stages respect provider coverage. This is field topology, not a particle
/// trajectory or plasma transport model. Refine the step and check convergence.
pub fn trace_field_line(
    field: &dyn MagneticField,
    start_m: Vec3,
    epoch: Epoch,
    step_m: f64,
    steps: usize,
) -> Result<Vec<Vec3>> {
    start_m.validate()?;
    finite(step_m, "field-line step")?;
    if step_m == 0.0 || steps == 0 {
        return Err(Error::InvalidInput(
            "field-line step and count must be nonzero".into(),
        ));
    }
    let direction = |point| -> Result<Vec3> { field.teslas_at(point, epoch)?.normalized() };
    let mut points = vec![start_m];
    let mut p = start_m;
    for _ in 0..steps {
        let k1 = direction(p)?;
        let k2 = direction(p + k1 * (step_m / 2.0))?;
        let k3 = direction(p + k2 * (step_m / 2.0))?;
        let k4 = direction(p + k3 * step_m)?;
        p = p + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (step_m / 6.0);
        p.validate()?;
        direction(p)?;
        points.push(p);
    }
    Ok(points)
}

#[derive(Debug, Clone, Copy)]
pub struct PlasmaIntegral {
    /// Integral ne ds, electrons/m²; 1 TECU = 1e16 electrons/m².
    pub electron_column_per_m2: f64,
    /// Integral ne B dot ds, signed along start -> end, SI units.
    pub electron_magnetic_column: f64,
}
impl PlasmaIntegral {
    /// First-order cold-plasma polarization rotation along the directed ray.
    /// Positive sign uses B dot ds; the receiver polarization convention may invert it.
    /// Assumes frequency well above plasma/cyclotron frequencies; not an HF ray tracer.
    pub fn faraday_rotation_rad(self, frequency_hz: f64) -> Result<f64> {
        finite(frequency_hz, "frequency")?;
        if frequency_hz <= 0.0 {
            return Err(Error::InvalidInput("frequency must be positive".into()));
        }
        finite(self.electron_magnetic_column, "magnetic column")?;
        let rotation = 2.3648e4 * self.electron_magnetic_column / frequency_hz / frequency_hz;
        finite(rotation, "Faraday rotation")?;
        Ok(rotation)
    }
}
/// Trapezoidal quadrature along a supplied straight segment and linearly
/// interpolated coordinate times. For a bent ray, sum separately integrated legs.
/// Neither this path nor its end times are inferred from the sampled plasma.
pub fn integrate_plasma(
    density: &dyn ElectronDensity,
    magnetic: &dyn MagneticField,
    start_m: Vec3,
    end_m: Vec3,
    start_epoch: Epoch,
    end_epoch: Epoch,
    segments: usize,
) -> Result<PlasmaIntegral> {
    start_m.validate()?;
    end_m.validate()?;
    if segments == 0 || end_epoch < start_epoch {
        return Err(Error::InvalidInput(
            "positive segment count and ordered epochs required".into(),
        ));
    }
    let displacement = end_m - start_m;
    let length = displacement.norm();
    finite(length, "path length")?;
    if length == 0.0 {
        return Err(Error::InvalidInput(
            "plasma path must have positive length".into(),
        ));
    }
    let direction = displacement * (1.0 / length);
    let dt = end_epoch.duration_since(start_epoch);
    let ds = length / segments as f64;
    let mut tec = 0.0;
    let mut magnetic_column = 0.0;
    for i in 0..=segments {
        let fraction = i as f64 / segments as f64;
        let point = start_m + displacement * fraction;
        let epoch = start_epoch.shifted(dt * fraction)?;
        let ne = density.per_m3_at(point, epoch)?;
        finite(ne, "electron density")?;
        if ne < 0.0 {
            return Err(Error::InvalidInput(
                "electron density must be nonnegative".into(),
            ));
        }
        let b = magnetic.teslas_at(point, epoch)?;
        b.validate()?;
        let weight = if i == 0 || i == segments { 0.5 } else { 1.0 };
        tec += ne * ds * weight;
        magnetic_column += ne * b.dot(direction) * ds * weight;
    }
    finite(tec, "electron column")?;
    finite(magnetic_column, "electron magnetic column")?;
    Ok(PlasmaIntegral {
        electron_column_per_m2: tec,
        electron_magnetic_column: magnetic_column,
    })
}
