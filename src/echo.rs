//! Coherent multi-target RF echoes at ONE common reception event.
//! Independent asteroids are not phase-locked transmitters. Their complex scattering
//! responses must be calibrated/modeled; unknown speckle is not recovered by geometry.
use crate::frames::{Epoch, State};
use crate::fusion::emission_epoch;
use crate::spacetime::Worldline;
use crate::{Error, Result, finite};
use std::collections::HashSet;

#[derive(Debug, Clone, Copy)]
pub struct ReceivedBistaticPath {
    pub emission_epoch: Epoch,
    pub scattering_epoch: Epoch,
    pub reception_epoch: Epoch,
    pub transmitter: State,
    pub target: State,
    pub receiver: State,
}
impl ReceivedBistaticPath {
    pub fn delay_s(self) -> f64 {
        self.reception_epoch.duration_since(self.emission_epoch)
    }
}
/// Both legs solved backward from a fixed reception epoch, unlike bistatic_path's
/// fixed TRANSMISSION epoch. Each scatterer may have a different emission epoch.
pub fn received_bistatic_path(
    transmitter: &dyn Worldline,
    target: &dyn Worldline,
    receiver: &dyn Worldline,
    reception_epoch: Epoch,
    tolerance_s: f64,
) -> Result<ReceivedBistaticPath> {
    let rx = receiver.state_at(reception_epoch)?;
    rx.position_m.validate()?;
    rx.velocity_m_s.validate()?;
    if rx.velocity_m_s.norm() >= crate::interferometry::SPEED_OF_LIGHT_M_S {
        return Err(Error::InvalidInput("receiver must be subluminal".into()));
    }
    let scattering_epoch = emission_epoch(target, rx.position_m, reception_epoch, tolerance_s)?;
    let scatter = target.state_at(scattering_epoch)?;
    let emitted = emission_epoch(
        transmitter,
        scatter.position_m,
        scattering_epoch,
        tolerance_s,
    )?;
    let tx = transmitter.state_at(emitted)?;
    tx.position_m.validate()?;
    tx.velocity_m_s.validate()?;
    scatter.position_m.validate()?;
    scatter.velocity_m_s.validate()?;
    if tx.velocity_m_s.norm() >= crate::interferometry::SPEED_OF_LIGHT_M_S
        || scatter.velocity_m_s.norm() >= crate::interferometry::SPEED_OF_LIGHT_M_S
    {
        return Err(Error::InvalidInput("path states must be subluminal".into()));
    }
    Ok(ReceivedBistaticPath {
        emission_epoch: emitted,
        scattering_epoch,
        reception_epoch,
        transmitter: tx,
        target: scatter,
        receiver: rx,
    })
}

#[derive(Debug, Clone, Copy, Default)]
pub struct ComplexAmplitude {
    pub real: f64,
    pub imaginary: f64,
}
impl ComplexAmplitude {
    pub fn rotated(self, phase_rad: f64) -> Result<Self> {
        finite(self.real, "complex amplitude")?;
        finite(self.imaginary, "complex amplitude")?;
        finite(phase_rad, "phase")?;
        let (sin, cos) = phase_rad.sin_cos();
        let result = Self {
            real: self.real * cos - self.imaginary * sin,
            imaginary: self.real * sin + self.imaginary * cos,
        };
        finite(result.real, "rotated amplitude")?;
        finite(result.imaginary, "rotated amplitude")?;
        Ok(result)
    }
    /// Squared amplitude, in the response amplitude's units squared (not implicitly watts).
    pub fn power(self) -> Result<f64> {
        let power = self.real * self.real + self.imaginary * self.imaginary;
        finite(power, "complex power")?;
        Ok(power)
    }
}

/// Response includes propagation amplitude, target spin/polarization, occlusion,
/// medium/instrument phase and waveform envelope as needed. Return zero for no echo.
/// Vacuum geometric carrier phase is applied separately. No implicit radar equation.
pub trait ScatterResponse: Send + Sync {
    fn amplitude(&self, path: &ReceivedBistaticPath, frequency_hz: f64)
    -> Result<ComplexAmplitude>;
}
impl<F> ScatterResponse for F
where
    F: Fn(&ReceivedBistaticPath, f64) -> Result<ComplexAmplitude> + Send + Sync,
{
    fn amplitude(
        &self,
        path: &ReceivedBistaticPath,
        frequency_hz: f64,
    ) -> Result<ComplexAmplitude> {
        self(path, frequency_hz)
    }
}
pub struct Scatterer<'a> {
    pub id: &'a str,
    pub trajectory: &'a dyn Worldline,
    pub response: &'a dyn ScatterResponse,
}
pub struct EchoContribution {
    pub id: String,
    pub path: ReceivedBistaticPath,
    pub amplitude: ComplexAmplitude,
}
pub struct CoherentEcho {
    pub channel: String,
    pub amplitude: ComplexAmplitude,
    pub contributions: Vec<EchoContribution>,
}
pub struct CwReference<'a> {
    pub channel: &'a str,
    pub frequency_hz: f64,
    pub emission_epoch: Epoch,
}
/// Sum calibrated echoes in one receiver/transmitter phase-linked CW channel.
/// Carrier phase is +2 pi f (emission - reference_emission); longer paths lag.
/// A nearby reference and split Epoch avoid multiplying carrier by astronomical dates.
/// Long light-times/high carriers still require a validated precision budget or
/// a specialized differential/relativistic solver. Sum unrelated channels as
/// likelihoods/powers, not by calling this function with arbitrary phase labels.
pub fn coherent_echo(
    reference: CwReference<'_>,
    transmitter: &dyn Worldline,
    receiver: &dyn Worldline,
    reception_epoch: Epoch,
    scatterers: &[Scatterer<'_>],
    tolerance_s: f64,
) -> Result<CoherentEcho> {
    let channel = reference.channel;
    let frequency_hz = reference.frequency_hz;
    let reference_emission = reference.emission_epoch;
    finite(frequency_hz, "carrier frequency")?;
    if channel.trim().is_empty() || frequency_hz <= 0.0 || scatterers.is_empty() {
        return Err(Error::InvalidInput(
            "channel, positive frequency and scatterers required".into(),
        ));
    }
    let mut ids = HashSet::new();
    let mut amplitude = ComplexAmplitude::default();
    let mut contributions = Vec::with_capacity(scatterers.len());
    for scatterer in scatterers {
        if scatterer.id.trim().is_empty() || !ids.insert(scatterer.id) {
            return Err(Error::InvalidInput(
                "scatterer IDs must be nonempty and unique".into(),
            ));
        }
        let path = received_bistatic_path(
            transmitter,
            scatterer.trajectory,
            receiver,
            reception_epoch,
            tolerance_s,
        )?;
        let cycles = frequency_hz * path.emission_epoch.duration_since(reference_emission);
        finite(cycles, "carrier cycles")?;
        let phase = std::f64::consts::TAU * cycles.rem_euclid(1.0);
        let echo = scatterer
            .response
            .amplitude(&path, frequency_hz)?
            .rotated(phase)?;
        amplitude.real += echo.real;
        amplitude.imaginary += echo.imaginary;
        contributions.push(EchoContribution {
            id: scatterer.id.into(),
            path,
            amplitude: echo,
        });
    }
    finite(amplitude.real, "summed echo")?;
    finite(amplitude.imaginary, "summed echo")?;
    Ok(CoherentEcho {
        channel: channel.into(),
        amplitude,
        contributions,
    })
}
