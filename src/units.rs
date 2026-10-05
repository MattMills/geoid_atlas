//! Explicit ingestion units. Core geometry remains SI; angles use named boundaries.
//! Quantities of different dimensions have distinct types. Variances must be
//! multiplied by the SQUARE of the unit factor, and cross covariances by both factors.
use crate::{Result, finite};

macro_rules! quantity {
    ($quantity:ident, $unit:ident, $getter:ident, {$($variant:ident => $factor:expr),+ $(,)?}) => {
        #[derive(Debug, Clone, Copy, PartialEq)]
        pub enum $unit { $($variant),+ }
        impl $unit {
            pub fn si_factor(self) -> f64 { match self { $(Self::$variant => $factor),+ } }
        }
        #[derive(Debug, Clone, Copy, PartialEq)]
        pub struct $quantity(f64);
        impl $quantity {
            pub fn new(value: f64, unit: $unit) -> Result<Self> {
                finite(value, "quantity")?;
                let si = value * unit.si_factor();
                finite(si, "SI quantity")?;
                Ok(Self(si))
            }
            pub fn $getter(self) -> f64 { self.0 }
            pub fn in_unit(self, unit: $unit) -> Result<f64> {
                let value = self.0 / unit.si_factor();
                finite(value, "converted quantity")?;
                Ok(value)
            }
        }
    }
}
quantity!(Length, LengthUnit, metres, {
    Metre => 1.0, Kilometre => 1_000.0, Centimetre => 0.01,
    InternationalFoot => 0.3048, UsSurveyFoot => 1200.0 / 3937.0,
    AstronomicalUnit => 149_597_870_700.0,
    Parsec => 149_597_870_700.0 * 648_000.0 / std::f64::consts::PI,
});
quantity!(Angle, AngleUnit, radians, {
    Radian => 1.0, Degree => std::f64::consts::PI / 180.0,
    Arcsecond => std::f64::consts::PI / 648_000.0,
    HourAngle => std::f64::consts::PI / 12.0,
});
quantity!(Duration, DurationUnit, seconds, {
    Second => 1.0, Millisecond => 1e-3, Microsecond => 1e-6,
    Nanosecond => 1e-9, Day => 86_400.0, JulianYear => 31_557_600.0,
});
quantity!(Frequency, FrequencyUnit, hertz, {
    Hertz => 1.0, Kilohertz => 1e3, Megahertz => 1e6, Gigahertz => 1e9,
});
quantity!(Pressure, PressureUnit, pascals, {
    Pascal => 1.0, Hectopascal => 100.0, StandardAtmosphere => 101_325.0,
});
quantity!(MagneticFluxDensity, MagneticUnit, teslas, {
    Tesla => 1.0, Microtesla => 1e-6, Nanotesla => 1e-9, Gauss => 1e-4,
});

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TemperatureUnit {
    Kelvin,
    Celsius,
    Fahrenheit,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Temperature(f64);
impl Temperature {
    pub fn new(value: f64, unit: TemperatureUnit) -> Result<Self> {
        finite(value, "temperature")?;
        let kelvin = match unit {
            TemperatureUnit::Kelvin => value,
            TemperatureUnit::Celsius => value + 273.15,
            TemperatureUnit::Fahrenheit => (value - 32.0) * (5.0 / 9.0) + 273.15,
        };
        finite(kelvin, "temperature")?;
        if kelvin < 0.0 {
            return Err(crate::Error::InvalidInput(
                "temperature below absolute zero".into(),
            ));
        }
        Ok(Self(kelvin))
    }
    pub fn kelvin(self) -> f64 {
        self.0
    }
    pub fn in_unit(self, unit: TemperatureUnit) -> Result<f64> {
        let value = match unit {
            TemperatureUnit::Kelvin => self.0,
            TemperatureUnit::Celsius => self.0 - 273.15,
            TemperatureUnit::Fahrenheit => (self.0 - 273.15) * (9.0 / 5.0) + 32.0,
        };
        finite(value, "converted temperature")?;
        Ok(value)
    }
}
