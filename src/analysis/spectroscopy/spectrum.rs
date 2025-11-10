use pyo3::prelude::*;
use crate::analysis::spectroscopy::tools::{hann_window, fft, alpha_w, mu_w, q_w};

#[pyfunction]
pub fn compute_spectra(
    mean: f64,
    vvacf: Vec<f64>,
    nfreq: usize,
    tcfl: usize,
    dt: f64,
    alpha_0: f64,
    mu_0: f64,
    temp: f64,
    spectrum: Option<String> // "IR" or "SFG"
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut vdos: Vec<f64> = Vec::new();

    let spec_type = spectrum.unwrap_or_else(|| "IR".to_string());
    let is_sfg = spec_type.eq_ignore_ascii_case("SFG");
    match spec_type.as_str() {
        "SFG" => println!("→ Calculating SFG spectrum"),
        "IR" => println!("→ Calculating IR spectrum"),
        _ => panic!("spectrum must be 'IR' or 'SFG'"),
    }

    // Apply the Hann window to the VVACF
    let vvacf_hann = hann_window(tcfl, dt, mean, &vvacf);

    

    if is_sfg{
        // Perform the FFT
        let (mut fft_img, mut fft_real) = fft(nfreq, tcfl, dt, &vvacf_hann);
        fft_img.iter_mut().for_each(|x| *x *= mean);
        fft_real.iter_mut().for_each(|x| *x *= mean);

        // Apply the quantum correction
        fft_img.iter_mut().enumerate().skip(1).for_each(|(freq, x)|
            *x *= q_w(freq as f64, temp)/(freq as f64).powi(2));
        fft_real.iter_mut().enumerate().skip(1).for_each(|(freq, x)|
            *x *= q_w(freq as f64, temp)/(freq as f64).powi(2));

        // Apply the non-condon effect
        fft_img.iter_mut().enumerate().skip(1).for_each(|(freq, x)|
            *x *= (alpha_w(freq as f64, alpha_0) * mu_w(freq as f64, mu_0)));
        fft_real.iter_mut().enumerate().skip(1).for_each(|(freq, x)|
            *x *= (alpha_w(freq as f64, alpha_0) * mu_w(freq as f64, mu_0)));

        (vvacf_hann, vdos, fft_real, fft_img)

    } else{
        let (mut fft_real, mut fft_img) = fft(nfreq, tcfl, dt, &vvacf_hann);

        vdos = fft_real.clone();
        // Apply the quantum correction
        fft_img.iter_mut().enumerate().skip(1).for_each(|(freq, x)|
            *x *= q_w(freq as f64, temp)/(freq as f64).powi(2));
        fft_real.iter_mut().enumerate().skip(1).for_each(|(freq, x)|
            *x *= q_w(freq as f64, temp)/(freq as f64).powi(2));

        // Apply the non-condon effect
        fft_img.iter_mut().enumerate().skip(1).for_each(|(freq, x)|
            *x *= mu_w(freq as f64, mu_0).powi(2));
        fft_real.iter_mut().enumerate().skip(1).for_each(|(freq, x)|
            *x *= mu_w(freq as f64, mu_0).powi(2));

        (vvacf_hann, vdos, fft_real, fft_img)
    } 

    // (vvacf_hann, vdos, fft_real, fft_img)

}