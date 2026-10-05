//! Uniform astronomical time coordinates with split Julian dates.
//!
//! UTC is deliberately excluded: use an authoritative leap-second table to
//! convert UTC before entering the model. TT = TAI + 32.184 s is supported;
//! TDB/TCB conversions require external relativistic ephemeris models.
use crate::frames::Epoch;
use crate::{Error, Result, finite};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimeScale {
    Tai,
    Tt,
    Tdb,
    Tcb,
}

/// Integral Julian day plus seconds from its noon boundary. Splitting the date
/// avoids losing sub-millisecond precision by subtracting large Julian floats.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JulianDate {
    day: i64,
    seconds: f64,
}
impl JulianDate {
    pub fn new(day: i64, seconds: f64) -> Result<Self> {
        finite(seconds, "Julian seconds")?;
        let shift = (seconds / 86_400.0).floor();
        if shift.abs() > 1e12 {
            return Err(Error::InvalidInput("Julian offset is too large".into()));
        }
        let day = day
            .checked_add(shift as i64)
            .ok_or_else(|| Error::InvalidInput("Julian day overflow".into()))?;
        let normalized = seconds.rem_euclid(86_400.0);
        if normalized >= 86_400.0 {
            return Ok(Self {
                day: day
                    .checked_add(1)
                    .ok_or_else(|| Error::InvalidInput("Julian day overflow".into()))?,
                seconds: 0.0,
            });
        }
        Ok(Self {
            day,
            seconds: normalized,
        })
    }
    pub fn day(self) -> i64 {
        self.day
    }
    pub fn seconds(self) -> f64 {
        self.seconds
    }
    pub fn shifted(self, seconds: f64) -> Result<Self> {
        Self::new(self.day, self.seconds + seconds)
    }
    pub fn difference_s(self, other: Self) -> Result<f64> {
        let days = self
            .day
            .checked_sub(other.day)
            .ok_or_else(|| Error::InvalidInput("Julian difference overflow".into()))?;
        let seconds = days as f64 * 86_400.0 + self.seconds - other.seconds;
        finite(seconds, "time difference")?;
        Ok(seconds)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Timestamp {
    pub date: JulianDate,
    pub scale: TimeScale,
}
impl Timestamp {
    pub fn to_scale(self, scale: TimeScale) -> Result<Self> {
        let delta = match (self.scale, scale) {
            (a, b) if a == b => 0.0,
            (TimeScale::Tai, TimeScale::Tt) => 32.184,
            (TimeScale::Tt, TimeScale::Tai) => -32.184,
            _ => {
                return Err(Error::InvalidInput(
                    "TDB/TCB conversion requires an external time model".into(),
                ));
            }
        };
        Ok(Self {
            date: self.date.shifted(delta)?,
            scale,
        })
    }
}

/// Epoch seconds are measured from this origin in this scale.
#[derive(Debug, Clone, Copy)]
pub struct TimeReference {
    pub origin: JulianDate,
    pub scale: TimeScale,
}
impl TimeReference {
    pub fn epoch_at(self, timestamp: Timestamp) -> Result<Epoch> {
        let date = timestamp.to_scale(self.scale)?.date;
        let days = date
            .day
            .checked_sub(self.origin.day)
            .ok_or_else(|| Error::InvalidInput("Julian difference overflow".into()))?;
        let whole = days
            .checked_mul(86_400)
            .ok_or_else(|| Error::InvalidInput("epoch overflow".into()))?;
        Epoch::from_parts(whole, date.seconds - self.origin.seconds)
    }
    pub fn timestamp_at(self, epoch: Epoch) -> Result<Timestamp> {
        let days = epoch.whole_seconds().div_euclid(86_400);
        let day = self
            .origin
            .day
            .checked_add(days)
            .ok_or_else(|| Error::InvalidInput("Julian day overflow".into()))?;
        let seconds = self.origin.seconds
            + epoch.whole_seconds().rem_euclid(86_400) as f64
            + epoch.fractional_seconds();
        Ok(Timestamp {
            date: JulianDate::new(day, seconds)?,
            scale: self.scale,
        })
    }
}
