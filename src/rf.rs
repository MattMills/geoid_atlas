//! Radio propagation baselines. These are explicit hypotheses, not inferred observations.
//! Friis assumes far-field free-space propagation; medium attenuation is homogeneous.
use crate::interferometry::SPEED_OF_LIGHT_M_S;
use crate::{Error, Result, finite};

pub fn free_space_path_loss_db(distance_m: f64, frequency_hz: f64) -> Result<f64> {
    finite(distance_m, "distance")?;
    finite(frequency_hz, "frequency")?;
    if distance_m <= 0.0 || frequency_hz <= 0.0 {
        return Err(Error::InvalidInput(
            "distance and frequency must be positive".into(),
        ));
    }
    // Sum logarithms to avoid overflow in large astronomical distance * frequency.
    Ok(
        20.0 * (4.0 * std::f64::consts::PI / SPEED_OF_LIGHT_M_S).log10()
            + 20.0 * distance_m.log10()
            + 20.0 * frequency_hz.log10(),
    )
}

#[derive(Debug, Clone, Copy)]
pub struct LinkBudget {
    pub transmitted_power_dbm: f64,
    pub transmit_gain_dbi: f64,
    pub receive_gain_dbi: f64,
    pub additional_loss_db: f64,
    pub shadowing_std_db: f64,
}
impl LinkBudget {
    pub fn predict(self, distance_m: f64, frequency_hz: f64) -> Result<PowerDistribution> {
        for v in [
            self.transmitted_power_dbm,
            self.transmit_gain_dbi,
            self.receive_gain_dbi,
            self.additional_loss_db,
        ] {
            finite(v, "link parameter")?;
        }
        if self.additional_loss_db < 0.0 {
            return Err(Error::InvalidInput(
                "additional loss must be nonnegative".into(),
            ));
        }
        PowerDistribution::new(
            self.transmitted_power_dbm + self.transmit_gain_dbi + self.receive_gain_dbi
                - self.additional_loss_db
                - free_space_path_loss_db(distance_m, frequency_hz)?,
            self.shadowing_std_db,
        )
    }
}

/// Gaussian received dBm, hence lognormal received watts. The median in watts
/// differs from its expectation when uncertainty is nonzero.
#[derive(Debug, Clone, Copy)]
pub struct PowerDistribution {
    mean_dbm: f64,
    std_db: f64,
}
#[derive(Debug, Clone, Copy, Default)]
pub struct PowerMoments {
    pub mean_w: f64,
    pub variance_w2: f64,
}
impl PowerDistribution {
    pub fn new(mean_dbm: f64, std_db: f64) -> Result<Self> {
        finite(mean_dbm, "mean power")?;
        finite(std_db, "power uncertainty")?;
        if std_db < 0.0 {
            return Err(Error::InvalidInput(
                "uncertainty must be nonnegative".into(),
            ));
        }
        let value = Self { mean_dbm, std_db };
        value.moments()?;
        Ok(value)
    }
    pub fn mean_dbm(self) -> f64 {
        self.mean_dbm
    }
    pub fn std_db(self) -> f64 {
        self.std_db
    }
    pub fn median_w(self) -> f64 {
        10.0_f64.powf((self.mean_dbm - 30.0) / 10.0)
    }
    pub fn moments(self) -> Result<PowerMoments> {
        let k = std::f64::consts::LN_10 / 10.0;
        let mu = k * (self.mean_dbm - 30.0);
        let variance = (k * self.std_db).powi(2);
        let mean = (mu + variance / 2.0).exp();
        let v = if variance == 0.0 {
            0.0
        } else {
            (2.0 * mu + variance).exp() * variance.exp_m1()
        };
        finite(mean, "mean watts")?;
        finite(v, "power variance")?;
        if mean == 0.0 {
            return Err(Error::InvalidInput(
                "power underflows representable range".into(),
            ));
        }
        Ok(PowerMoments {
            mean_w: mean,
            variance_w2: v,
        })
    }
}
/// Incoherent sum of mutually independent signals, not coherent field interference.
/// Caller must establish independence and frequency-band compatibility.
pub fn independent_power_sum(signals: &[PowerDistribution]) -> Result<PowerMoments> {
    let mut total = PowerMoments::default();
    for s in signals {
        let m = s.moments()?;
        total.mean_w += m.mean_w;
        total.variance_w2 += m.variance_w2;
    }
    finite(total.mean_w, "total mean power")?;
    finite(total.variance_w2, "total variance")?;
    Ok(total)
}

/// Homogeneous isotropic material with real relative permittivity/permeability,
/// and Ohmic conductivity. Frequency-dependent material data must be supplied
/// independently for each evaluation; interface refraction and scattering are absent.
pub struct Material {
    pub relative_permittivity: f64,
    pub relative_permeability: f64,
    pub conductivity_s_m: f64,
}
#[derive(Debug, Clone, Copy)]
pub struct PropagationConstant {
    pub attenuation_nepers_m: f64,
    pub phase_rad_m: f64,
}
impl Material {
    pub fn propagation_constant(&self, frequency_hz: f64) -> Result<PropagationConstant> {
        for v in [
            frequency_hz,
            self.relative_permittivity,
            self.relative_permeability,
            self.conductivity_s_m,
        ] {
            finite(v, "material parameter")?;
        }
        if frequency_hz <= 0.0
            || self.relative_permittivity < 1.0
            || self.relative_permeability <= 0.0
            || self.conductivity_s_m < 0.0
        {
            return Err(Error::InvalidInput(
                "require f>0, relative permittivity>=1, permeability>0, conductivity>=0".into(),
            ));
        }
        let omega = 2.0 * std::f64::consts::PI * frequency_hz;
        let mu0 = 1.256_637_061_27e-6;
        let eps0 = 1.0 / (mu0 * SPEED_OF_LIGHT_M_S.powi(2));
        let eps = eps0 * self.relative_permittivity;
        let mu = mu0 * self.relative_permeability;
        let ratio = self.conductivity_s_m / (omega * eps);
        let root = ratio.hypot(1.0);
        let base = omega * (mu * eps / 2.0).sqrt();
        let alpha = base * ratio / (root + 1.0).sqrt();
        let beta = base * (root + 1.0).sqrt();
        finite(alpha, "attenuation")?;
        finite(beta, "phase constant")?;
        Ok(PropagationConstant {
            attenuation_nepers_m: alpha,
            phase_rad_m: beta,
        })
    }
}
impl PropagationConstant {
    pub fn absorption_loss_db(self, length_m: f64) -> Result<f64> {
        finite(length_m, "path length")?;
        finite(self.attenuation_nepers_m, "attenuation")?;
        if length_m < 0.0 || self.attenuation_nepers_m < 0.0 {
            return Err(Error::InvalidInput(
                "path and attenuation must be nonnegative".into(),
            ));
        }
        let loss = 20.0 / std::f64::consts::LN_10 * self.attenuation_nepers_m * length_m;
        finite(loss, "absorption loss")?;
        Ok(loss)
    }
}

/// Approximate radio tropospheric refractive index from temperature and total
/// and water-vapour partial pressures. Pressures are Pa, temperature is K.
/// Uses `N = 77.6 p[hPa]/T + 3.73e5 e[hPa]/T²`. Not an optical dispersion law.
pub fn tropospheric_refractive_index(
    temperature_k: f64,
    total_pressure_pa: f64,
    water_vapour_pressure_pa: f64,
) -> Result<f64> {
    for v in [temperature_k, total_pressure_pa, water_vapour_pressure_pa] {
        finite(v, "weather parameter")?;
    }
    if temperature_k <= 0.0
        || total_pressure_pa <= 0.0
        || water_vapour_pressure_pa < 0.0
        || water_vapour_pressure_pa > total_pressure_pa
    {
        return Err(Error::InvalidInput(
            "invalid temperature or partial pressures".into(),
        ));
    }
    let refractivity = 77.6 * (total_pressure_pa / 100.0) / temperature_k
        + 3.73e5 * (water_vapour_pressure_pa / 100.0) / temperature_k.powi(2);
    let index = 1.0 + refractivity * 1e-6;
    finite(index, "refractive index")?;
    Ok(index)
}

#[derive(Debug, Clone, Copy)]
pub struct IonosphericDelay {
    pub group_delay_s: f64,
    pub phase_delay_s: f64,
}
/// First-order cold-plasma correction at frequencies well above plasma cutoff.
/// TEC is the ray-integrated electron column in electrons/m² (1 TECU = 1e16).
/// No magnetic-field, Faraday rotation, or higher-order plasma terms are implicit.
pub fn ionospheric_delay(
    total_electron_content_per_m2: f64,
    frequency_hz: f64,
) -> Result<IonosphericDelay> {
    finite(total_electron_content_per_m2, "TEC")?;
    finite(frequency_hz, "frequency")?;
    if total_electron_content_per_m2 < 0.0 || frequency_hz <= 0.0 {
        return Err(Error::InvalidInput("require TEC>=0 and f>0".into()));
    }
    let delay =
        40.3 / SPEED_OF_LIGHT_M_S * total_electron_content_per_m2 / frequency_hz / frequency_hz;
    finite(delay, "ionospheric delay")?;
    Ok(IonosphericDelay {
        group_delay_s: delay,
        phase_delay_s: -delay,
    })
}
