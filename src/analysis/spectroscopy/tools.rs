use rayon::prelude::*;
use std::f64::consts::PI;
use crate::utils::physchem::constants::{C, H_BAR, K_B};



pub fn q_w(freq: f64, temp: f64) -> f64 {
    let c = C / 1.0e-15; // Adjusted speed of light
    let beta = 1.0 / (K_B * temp); // Inverse thermal energy
    let exponent = -beta * H_BAR * freq * c;
    let q_w = beta * H_BAR * freq * c / (1.0 - exponent.exp()); // Compute Q_w
    q_w
}

pub fn mu_w(freq: f64, mu_0: f64) -> f64 {
    let mut mu_w = 1.377 + (3737.0 - freq) * 53.03 / 6932.2;
    mu_w *= mu_0; // Multiply by Mu_0
    mu_w
}

pub fn alpha_w(freq: f64, alpha_0: f64) -> f64 {
    let mut alpha_w = 1.271 + (3737.0 - freq) * 6.287 / 6932.2;
    alpha_w *= alpha_0; // Multiply by Alpha_0
    alpha_w
}

pub fn norm(position:[f64;3]) -> f64{
    let dist = (position[0].powi(2) + position[1].powi(2) + position[2].powi(2)).sqrt();
    dist
}

pub fn dp(position:[f64;3], velocity:[f64;3]) -> f64{
    let projection = position[0]*velocity[0] + position[1]*velocity[1] + position[2]*velocity[2];
    projection
}

pub fn hann_window(tcfl:usize, dt:f64, mean:f64, tcf: &[f64]) -> Vec<f64> {
    let tcf_hann = (0..tcfl)
                .into_par_iter()
                .map(|t_| {
                    let hann_factor = (PI * t_ as f64 * dt / (tcfl as f64 * dt * 2.0)).cos().powi(2);
                    tcf[t_] * hann_factor / mean
                })
                .collect::<Vec<f64>>();

    tcf_hann
}

pub fn fft(nfreq: usize, tcfl: usize, dt: f64, tcf_hann: &[f64]) -> (Vec<f64>, Vec<f64>){
    let fft_results: Vec<(f64, f64)> = (0..nfreq)
        .into_par_iter()
        .map(|freq| {
            let mut real_sum = 0.0;
            let mut imag_sum = 0.0;

            for t_ in 0..tcfl {
                let angle = 2.0 * PI * t_ as f64 * dt * freq as f64 * C;
                real_sum += tcf_hann[t_] * angle.cos() * dt;
                imag_sum += tcf_hann[t_] * angle.sin() * dt;
            }
            (real_sum, imag_sum)
        })
        .collect();

    let real: Vec<f64> = fft_results.iter().map(|&(re, _)| re).collect();
    let imag: Vec<f64> = fft_results.iter().map(|&(_, im)| im).collect();

    (real, imag)
}

