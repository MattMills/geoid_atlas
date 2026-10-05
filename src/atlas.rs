//! Body-scoped source catalogue and scalar expectation atlas.
//!
//! Observations and model hypotheses remain distinguishable. Fusion assumes
//! independent, unbiased Gaussian errors BETWEEN named independence groups;
//! within a group only the most precise available estimate is used.
use crate::coordinates::{Crs, Ellipsoid, Geodetic};
use crate::frames::{Body, Epoch, FrameGraph, State, Vec3};
use crate::gravity::{GravityField, PointMass, ellipsoidal_height, point_mass_field};
use crate::raster::{Drg, Raster, Rgb};
use crate::{Error, Result, finite};
use std::collections::{BTreeMap, HashMap};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceKind {
    Observation,
    Model,
}
#[derive(Debug, Clone)]
pub struct Source {
    pub id: String,
    pub title: String,
    /// Dataset URI, citation, or reproducible model specification. Not fetched implicitly.
    pub reference: String,
    pub kind: SourceKind,
    /// Shared error lineage, e.g. two products derived from the same survey.
    pub independence_group: String,
}
#[derive(Debug, Clone, Copy)]
pub enum EvidenceSelection {
    Observed,
    Modeled,
    All,
}
impl EvidenceSelection {
    fn includes(self, kind: SourceKind) -> bool {
        match self {
            Self::Observed => kind == SourceKind::Observation,
            Self::Modeled => kind == SourceKind::Model,
            Self::All => true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Band {
    lower_hz: f64,
    upper_hz: f64,
}
impl Band {
    pub fn new(lower_hz: f64, upper_hz: f64) -> Result<Self> {
        finite(lower_hz, "band lower bound")?;
        finite(upper_hz, "band upper bound")?;
        if lower_hz <= 0.0 || upper_hz <= lower_hz {
            return Err(Error::InvalidInput(
                "band requires 0 < lower < upper".into(),
            ));
        }
        Ok(Self { lower_hz, upper_hz })
    }
    pub fn lower_hz(self) -> f64 {
        self.lower_hz
    }
    pub fn upper_hz(self) -> f64 {
        self.upper_hz
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Quantity {
    EllipsoidalElevationM,
    OrthometricElevationM,
    GeoidUndulationM,
    GravityAccelerationMPerS2,
    RefractiveIndex,
    ConductivitySPerM,
    TemperatureK,
    PressurePa,
    RelativeHumidityFraction,
    ElectronDensityPerM3,
    /// Band-integrated received power in dBm; fuse only identical band definitions.
    ReceivedPowerDbm(Band),
}

#[derive(Debug, Clone, Copy)]
pub struct Estimate {
    value: f64,
    std_dev: f64,
}
impl Estimate {
    pub fn new(value: f64, std_dev: f64) -> Result<Self> {
        finite(value, "estimate")?;
        finite(std_dev, "standard deviation")?;
        if std_dev <= 0.0 {
            return Err(Error::InvalidInput(
                "standard deviation must be positive".into(),
            ));
        }
        Ok(Self { value, std_dev })
    }
    pub fn value(self) -> f64 {
        self.value
    }
    pub fn std_dev(self) -> f64 {
        self.std_dev
    }
}

/// Extension point for tiles, files, measured samples, material volumes, or external models.
/// Geographic coordinates and elevations are relative to the associated body's ellipsoid.
pub trait ScalarField: Send + Sync {
    fn sample(&self, point: Geodetic, epoch: Epoch) -> Result<Estimate>;
}
pub struct ConstantField {
    pub estimate: Estimate,
}
impl ScalarField for ConstantField {
    fn sample(&self, _: Geodetic, _: Epoch) -> Result<Estimate> {
        Ok(self.estimate)
    }
}
pub struct RasterField {
    pub raster: Raster<f64>,
    pub std_dev: f64,
}
impl ScalarField for RasterField {
    fn sample(&self, point: Geodetic, _: Epoch) -> Result<Estimate> {
        Estimate::new(self.raster.bilinear(point)?, self.std_dev)
    }
}

/// Inclusive validity interval in the graph's declared coordinate time scale.
#[derive(Debug, Clone, Copy)]
pub struct Validity {
    start: Epoch,
    end: Epoch,
}
impl Validity {
    pub fn new(start: Epoch, end: Epoch) -> Result<Self> {
        if start > end {
            return Err(Error::InvalidInput("validity start exceeds end".into()));
        }
        Ok(Self { start, end })
    }
    pub fn contains(self, epoch: Epoch) -> bool {
        self.start <= epoch && epoch <= self.end
    }
}
struct Layer {
    body: String,
    source: String,
    quantity: Quantity,
    validity: Option<Validity>,
    field: Box<dyn ScalarField>,
}
struct ImageLayer {
    body: String,
    source: String,
    validity: Option<Validity>,
    image: Drg,
}

#[derive(Debug, Clone)]
pub struct AtlasEstimate {
    pub estimate: Estimate,
    pub contributing_sources: Vec<String>,
    pub redundant_sources: Vec<String>,
    /// Diagnostic inconsistency under the independent Gaussian assumption.
    pub reduced_chi_squared: Option<f64>,
}

pub struct Atlas {
    frames: FrameGraph,
    bodies: BTreeMap<String, Body>,
    sources: HashMap<String, Source>,
    layers: Vec<Layer>,
    images: Vec<ImageLayer>,
}
impl Atlas {
    pub fn new(frames: FrameGraph) -> Self {
        Self {
            frames,
            bodies: BTreeMap::new(),
            sources: HashMap::new(),
            layers: vec![],
            images: vec![],
        }
    }
    pub fn frames(&self) -> &FrameGraph {
        &self.frames
    }
    pub fn add_body(&mut self, id: impl Into<String>, body: Body) -> Result<()> {
        let id = id.into();
        self.frames.name(body.fixed_frame)?;
        finite(body.gravitational_parameter_m3_s2, "body mu")?;
        if id.trim().is_empty()
            || body.name.trim().is_empty()
            || self.bodies.contains_key(&id)
            || body.gravitational_parameter_m3_s2 < 0.0
        {
            return Err(Error::InvalidInput(
                "body needs unique ID, name, and nonnegative mu".into(),
            ));
        }
        self.bodies.insert(id, body);
        Ok(())
    }
    pub fn body(&self, id: &str) -> Result<&Body> {
        self.bodies
            .get(id)
            .ok_or_else(|| Error::InvalidInput(format!("unknown body: {id}")))
    }
    pub fn add_source(&mut self, source: Source) -> Result<()> {
        if source.id.trim().is_empty()
            || source.title.trim().is_empty()
            || source.reference.trim().is_empty()
            || source.independence_group.trim().is_empty()
        {
            return Err(Error::InvalidInput(
                "source metadata must be nonempty".into(),
            ));
        }
        if self.sources.contains_key(&source.id) {
            return Err(Error::DuplicateSource(source.id));
        }
        self.sources.insert(source.id.clone(), source);
        Ok(())
    }
    pub fn source(&self, id: &str) -> Result<&Source> {
        self.sources
            .get(id)
            .ok_or_else(|| Error::UnknownSource(id.into()))
    }
    pub fn add_scalar(
        &mut self,
        body: &str,
        source: &str,
        quantity: Quantity,
        validity: Option<Validity>,
        field: impl ScalarField + 'static,
    ) -> Result<()> {
        self.body(body)?;
        self.source(source)?;
        if self
            .layers
            .iter()
            .any(|l| l.body == body && l.source == source && l.quantity == quantity)
        {
            return Err(Error::InvalidInput(
                "one scalar field per body/source/quantity; tile within a field adapter".into(),
            ));
        }
        self.layers.push(Layer {
            body: body.into(),
            source: source.into(),
            quantity,
            validity,
            field: Box::new(field),
        });
        Ok(())
    }
    /// Validates the Earth-only projection restriction before adding a raster field.
    pub fn add_raster(
        &mut self,
        body: &str,
        source: &str,
        quantity: Quantity,
        validity: Option<Validity>,
        field: RasterField,
    ) -> Result<()> {
        self.validate_crs(body, field.raster.crs())?;
        finite(field.std_dev, "raster uncertainty")?;
        if field.std_dev <= 0.0 {
            return Err(Error::InvalidInput(
                "raster uncertainty must be positive".into(),
            ));
        }
        self.add_scalar(body, source, quantity, validity, field)
    }
    fn validate_crs(&self, body: &str, crs: Crs) -> Result<()> {
        if self.body(body)?.ellipsoid != Ellipsoid::WGS84 && crs != Crs::Planetographic {
            return Err(Error::InvalidInput(
                "non-WGS84 bodies require Planetographic rasters or custom field reprojection"
                    .into(),
            ));
        }
        Ok(())
    }
    pub fn add_drg(
        &mut self,
        body: &str,
        source: &str,
        validity: Option<Validity>,
        image: Drg,
    ) -> Result<()> {
        self.source(source)?;
        self.validate_crs(body, image.crs())?;
        if self
            .images
            .iter()
            .any(|i| i.body == body && i.source == source)
        {
            return Err(Error::InvalidInput("duplicate body/source image".into()));
        }
        self.images.push(ImageLayer {
            body: body.into(),
            source: source.into(),
            validity,
            image,
        });
        Ok(())
    }
    pub fn drg_pixel(
        &self,
        body: &str,
        source: &str,
        point: Geodetic,
        epoch: Epoch,
    ) -> Result<Rgb> {
        self.body(body)?;
        self.source(source)?;
        let image = self
            .images
            .iter()
            .find(|i| {
                i.body == body && i.source == source && i.validity.is_none_or(|v| v.contains(epoch))
            })
            .ok_or(Error::NoData)?;
        Ok(*image.image.nearest(point)?)
    }
    pub fn expectation(
        &self,
        body: &str,
        quantity: Quantity,
        point: Geodetic,
        epoch: Epoch,
        selection: EvidenceSelection,
    ) -> Result<AtlasEstimate> {
        self.body(body)?;
        let mut groups: BTreeMap<&str, (&str, Estimate)> = BTreeMap::new();
        let mut redundant = vec![];
        for layer in &self.layers {
            let source = self.source(&layer.source)?;
            if layer.body != body
                || layer.quantity != quantity
                || !selection.includes(source.kind)
                || !layer.validity.is_none_or(|v| v.contains(epoch))
            {
                continue;
            }
            let estimate = match layer.field.sample(point, epoch) {
                Ok(e) => e,
                Err(Error::NoData | Error::OutsideCoverage) => continue,
                Err(e) => return Err(e),
            };
            match groups.get(source.independence_group.as_str()) {
                Some((_, old)) if old.std_dev <= estimate.std_dev => {
                    redundant.push(source.id.clone())
                }
                Some((old_id, _)) => {
                    redundant.push((*old_id).into());
                    groups.insert(&source.independence_group, (&source.id, estimate));
                }
                None => {
                    groups.insert(&source.independence_group, (&source.id, estimate));
                }
            }
        }
        if groups.is_empty() {
            return Err(Error::NoData);
        }
        let min_sigma = groups
            .values()
            .map(|(_, e)| e.std_dev)
            .fold(f64::INFINITY, f64::min);
        let total_weight: f64 = groups
            .values()
            .map(|(_, e)| (min_sigma / e.std_dev).powi(2))
            .sum();
        let mean: f64 = groups
            .values()
            .map(|(_, e)| e.value * ((min_sigma / e.std_dev).powi(2) / total_weight))
            .sum();
        let chi2: f64 = groups
            .values()
            .map(|(_, e)| ((e.value - mean) / e.std_dev).powi(2))
            .sum();
        finite(chi2, "evidence inconsistency")?;
        let reduced = if groups.len() > 1 {
            Some(chi2 / (groups.len() - 1) as f64)
        } else {
            None
        };
        // Inflate formal uncertainty for disagreement; report the diagnostic.
        let sigma = min_sigma / total_weight.sqrt() * reduced.unwrap_or(1.0).max(1.0).sqrt();
        redundant.sort();
        Ok(AtlasEstimate {
            estimate: Estimate::new(mean, sigma)?,
            contributing_sources: groups.values().map(|(id, _)| (*id).into()).collect(),
            redundant_sources: redundant,
            reduced_chi_squared: reduced,
        })
    }
    pub fn site_in_root(
        &self,
        body: &str,
        ellipsoidal_point: Geodetic,
        epoch: Epoch,
    ) -> Result<State> {
        let b = self.body(body)?;
        self.frames.transform(
            b.surface_state(ellipsoidal_point),
            b.fixed_frame,
            self.frames.root(),
            epoch,
        )
    }
    /// Orthometric DEM conversion uses explicitly selected elevation/geoid expectations.
    /// Returns the mean position only, not a propagated position covariance.
    pub fn terrain_site_in_root(
        &self,
        body: &str,
        point: Geodetic,
        epoch: Epoch,
        datum: crate::raster::VerticalDatum,
        selection: EvidenceSelection,
    ) -> Result<State> {
        let height = match datum {
            crate::raster::VerticalDatum::Ellipsoid => {
                self.expectation(
                    body,
                    Quantity::EllipsoidalElevationM,
                    point,
                    epoch,
                    selection,
                )?
                .estimate
                .value
            }
            crate::raster::VerticalDatum::Orthometric => ellipsoidal_height(
                self.expectation(
                    body,
                    Quantity::OrthometricElevationM,
                    point,
                    epoch,
                    selection,
                )?
                .estimate
                .value,
                self.expectation(body, Quantity::GeoidUndulationM, point, epoch, selection)?
                    .estimate
                    .value,
            )?,
        };
        self.site_in_root(body, point.with_height(height)?, epoch)
    }
    /// All registered masses evaluated at one epoch; unsuitable inside extended bodies.
    pub fn gravity_in_root(&self, point: Vec3, epoch: Epoch) -> Result<GravityField> {
        let masses = self
            .bodies
            .values()
            .map(|b| {
                Ok(PointMass {
                    position_m: self
                        .frames
                        .pose_in_root(b.fixed_frame, epoch)?
                        .origin
                        .position_m,
                    gravitational_parameter_m3_s2: b.gravitational_parameter_m3_s2,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        point_mass_field(point, &masses)
    }
}
