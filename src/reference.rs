//! Explicit vertical references and body-local tangent/centric coordinates.
//! Horizontal datum shifts need supplied transforms, not a changed ellipsoid label.
use crate::coordinates::{Ecef, Ellipsoid, Geodetic};
use crate::frames::{Epoch, FrameGraph, FrameId, Pose, Rotation, State, Vec3};
use crate::gravity::{GeoidModel, ellipsoidal_height, orthometric_height};
use crate::{Error, Result, finite};

pub enum HeightReference<'a> {
    Ellipsoidal,
    /// The caller supplies a model in this exact horizontal/body reference.
    Orthometric(&'a dyn GeoidModel),
}
pub struct GeodeticReference<'a> {
    pub frame: FrameId,
    pub ellipsoid: Ellipsoid,
    pub height: HeightReference<'a>,
}
/// Geographic input whose height follows the accompanying reference (not always h).
#[derive(Debug, Clone, Copy)]
pub struct GeographicCoordinate {
    pub latitude_deg: f64,
    pub longitude_deg: f64,
    pub height_m: f64,
}
impl GeodeticReference<'_> {
    pub fn to_cartesian(&self, coordinate: GeographicCoordinate) -> Result<Vec3> {
        let point = Geodetic::new(
            coordinate.latitude_deg,
            coordinate.longitude_deg,
            coordinate.height_m,
        )?;
        // Geoid evaluated on the reference surface; input H is not ellipsoidal height.
        let surface = point.with_height(0.0)?;
        let h = match self.height {
            HeightReference::Ellipsoidal => point.height_m(),
            HeightReference::Orthometric(model) => {
                ellipsoidal_height(point.height_m(), model.undulation_m(surface)?)?
            }
        };
        let p = self.ellipsoid.to_cartesian(point.with_height(h)?);
        Vec3::new(p.x, p.y, p.z)
    }
    pub fn from_cartesian(&self, position_m: Vec3) -> Result<GeographicCoordinate> {
        let point = self.ellipsoid.to_geodetic(Ecef {
            x: position_m.x,
            y: position_m.y,
            z: position_m.z,
        })?;
        let height_m = match self.height {
            HeightReference::Ellipsoidal => point.height_m(),
            HeightReference::Orthometric(model) => orthometric_height(
                point.height_m(),
                model.undulation_m(point.with_height(0.0)?)?,
            )?,
        };
        Ok(GeographicCoordinate {
            latitude_deg: point.latitude_deg(),
            longitude_deg: point.longitude_deg(),
            height_m,
        })
    }
    /// Represent the SAME event in another body's/datum's axes at one epoch.
    /// This does not move a station onto another body. The target model may lack coverage.
    pub fn convert(
        &self,
        coordinate: GeographicCoordinate,
        target: &Self,
        graph: &FrameGraph,
        epoch: Epoch,
    ) -> Result<GeographicCoordinate> {
        let state = graph.transform(
            State::stationary(self.to_cartesian(coordinate)?),
            self.frame,
            target.frame,
            epoch,
        )?;
        target.from_cartesian(state.position_m)
    }
}

/// ENU local-to-body pose on ANY oblate ellipsoid. At a pole, longitude chooses
/// the tangent axes; ENU longitude itself has no unique physical polar direction.
pub fn tangent_pose(ellipsoid: Ellipsoid, origin: Geodetic) -> Result<Pose> {
    let (sl, cl) = origin.latitude_deg().to_radians().sin_cos();
    let (so, co) = origin.longitude_deg().to_radians().sin_cos();
    let p = ellipsoid.to_cartesian(origin);
    Ok(Pose {
        origin: State::stationary(Vec3::new(p.x, p.y, p.z)?),
        rotation: Rotation::new([
            [-so, -sl * co, cl * co],
            [co, -sl * so, cl * so],
            [0.0, cl, sl],
        ])?,
        angular_velocity_rad_s: Vec3::ZERO,
    })
}

/// Body-centred east-positive longitude and geocentric latitude, radius in metres.
/// Planetographic latitude uses the ellipsoid normal and is generally different.
#[derive(Debug, Clone, Copy)]
pub struct Planetocentric {
    pub latitude_deg: f64,
    pub longitude_deg: f64,
    pub radius_m: f64,
}
impl Planetocentric {
    pub fn to_cartesian(self) -> Result<Vec3> {
        let angles = Geodetic::new(self.latitude_deg, self.longitude_deg, 0.0)?;
        finite(self.radius_m, "radius")?;
        if self.radius_m <= 0.0 {
            return Err(Error::InvalidInput("radius must be positive".into()));
        }
        let (sl, cl) = angles.latitude_deg().to_radians().sin_cos();
        let (so, co) = angles.longitude_deg().to_radians().sin_cos();
        Vec3::new(
            self.radius_m * cl * co,
            self.radius_m * cl * so,
            self.radius_m * sl,
        )
    }
    pub fn from_cartesian(p: Vec3) -> Result<Self> {
        p.validate()?;
        let radius_m = p.norm();
        finite(radius_m, "radius")?;
        if radius_m == 0.0 {
            return Err(Error::InvalidInput(
                "centric latitude undefined at centre".into(),
            ));
        }
        Ok(Self {
            latitude_deg: p.z.atan2(p.x.hypot(p.y)).to_degrees(),
            longitude_deg: p.y.atan2(p.x).to_degrees(),
            radius_m,
        })
    }
}
