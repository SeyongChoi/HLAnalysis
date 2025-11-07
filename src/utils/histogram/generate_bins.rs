use pyo3::prelude::*;
use crate::atoms::Atoms;


/// Generate 1D bin edges and their centers on [start, end).
///
/// # Arguments
/// * `start`    - Lower bound of the range (inclusive).
/// * `end`      - Upper bound of the range (exclusive).
/// * `interval` - Bin width.
///
/// # Returns
/// * `(bins, bins_centers)`
///   - `bins`:       Monotonic bin edges of length N (right-open intervals).
///   - `bins_centers`: Centers computed from adjacent edges (length N-1).
#[pyfunction]
pub fn generate_bins(start: f64, end: f64, interval: f64) -> (Vec<f64>, Vec<f64>) {
    let mut bins = Vec::new();
    let mut current = start;

    // Build right-open bins: [bins[i], bins[i+1])
    while current < end {
        bins.push(current);
        current += interval;
    }

    // Compute the centers of each bin: (left + right) / 2
    let bins_centers: Vec<f64> = bins.windows(2).map(|w| (w[0] + w[1]) / 2.0).collect();

    (bins, bins_centers)
}

/// Group per-molecule angles into z-bins along the surface normal (z-axis).
///
/// Each entry in `angles` is `(o_idx, h_idx, theta)`. We take the oxygen's z-position
/// from `atoms`, find which z-bin it falls into, then push `theta` to that bin.
///
/// # Arguments
/// * `atoms`          - Atom container providing positions.
/// * `angles`         - List of (oxygen index, hydrogen index, angle).
/// * `bins_normal`    - z bin edges (length K, right-open).
/// * `centers_normal` - z bin centers (length K-1).
///
/// # Returns
/// * `Vec<Vec<f64>>` of length `centers_normal.len()`, where each inner vector
///   stores all angles whose O atom falls into that z-bin.
pub fn generate_angle_bins_along_normal(
    atoms: &Atoms,
    angles: &Vec<(i32, i32, f64)>,
    bins_normal: &[f64],
    centers_normal: &[f64],
) -> Vec<Vec<f64>> {
    // Initialize a vector of empty angle lists, one per z-bin.
    let mut angle_bins: Vec<Vec<f64>> = vec![Vec::new(); centers_normal.len()];

    if let Some(positions) = &atoms.positions() {

        for &(o_idx, _h_idx, theta) in angles {
            let z_oxygen = positions[o_idx as usize][2];

            // Search for the right-open interval [bins[i], bins[i+1])
            let mut bin_idx = None;
            for i in 0..(bins_normal.len() - 1) {
                if bins_normal[i] <= z_oxygen && z_oxygen < bins_normal[i + 1] {
                    bin_idx = Some(i);
                    break;
                }
            }

            // Append the angle to the matched z-bin
            if let Some(idx) = bin_idx {
                angle_bins[idx].push(theta);
            }
        }
    }
    angle_bins
}


/// Same as `generate_angle_bins_along_normal` but semantically dedicated to OH angles.
///
/// This duplicates the structure for clarity with OH-specific pipelines,
/// while keeping identical logic.
///
/// # Arguments
/// * `atoms`          - Atom container providing positions.
/// * `angles`         - List of (oxygen index, hydrogen index, angle).
/// * `bins_normal`    - z bin edges (length K, right-open).
/// * `centers_normal` - z bin centers (length K-1).
///
/// # Returns
/// * `Vec<Vec<f64>>` of per-z-bin angle lists.
pub fn generate_oh_angle_bins_along_normal(
    atoms: &Atoms,
    angles: &Vec<(i32, i32, f64)>,
    bins_normal: &[f64],
    centers_normal: &[f64],
) -> Vec<Vec<f64>> {
    
    let mut angle_bins: Vec<Vec<f64>> = vec![Vec::new(); centers_normal.len()];

    if let Some(positions) = &atoms.positions() {
    
        for &(o_idx, _h_idx, theta) in angles {
            let z_oxygen = positions[o_idx as usize][2];

            let mut bin_idx = None;
            for i in 0..(bins_normal.len() - 1) {
                if bins_normal[i] <= z_oxygen && z_oxygen < bins_normal[i + 1] {
                    bin_idx = Some(i);
                    break;
                }
            }

            // If z_oxygen is within the range of bins, assign to the correct bin
            if let Some(idx) = bin_idx {
                angle_bins[idx].push(theta);
            }
        }
    }
    angle_bins
}

/// Bin dipole-related angles by the oxygen atom's z-position.
///
/// Each entry in `angles` is `(o_idx, h1_idx, h2_idx, theta)`. We read O's z,
/// locate its z-bin, then push `theta`.
///
/// # Arguments
/// * `atoms`          - Atom container providing positions.
/// * `angles`         - List of (oxygen index, H1 index, H2 index, angle).
/// * `bins_normal`    - z bin edges (length K, right-open).
/// * `centers_normal` - z bin centers (length K-1).
///
/// # Returns
/// * `Vec<Vec<f64>>` of per-z-bin angle lists.
pub fn generate_dipole_angle_bins_along_normal(
    atoms: &Atoms,
    angles: &Vec<(i32, i32, i32, f64)>,
    bins_normal: &[f64],
    centers_normal: &[f64],
) -> Vec<Vec<f64>> {

    let mut angle_bins: Vec<Vec<f64>> = vec![Vec::new(); centers_normal.len()];

    if let Some(positions) = &atoms.positions() {
        for &(o_idx, _h1_idx, _h2_idx, theta) in angles {
            let z_oxygen = positions[o_idx as usize][2];

            let mut bin_idx = None;
            for i in 0..(bins_normal.len() - 1) {
                if bins_normal[i] <= z_oxygen && z_oxygen < bins_normal[i + 1] {
                    bin_idx = Some(i);
                    break;
                }
            }

            if let Some(idx) = bin_idx {
                angle_bins[idx].push(theta);
            }
        }
    }
    angle_bins
}

/// Bin angle pairs (dipole angle, H–H angle) by the oxygen atom's z-position.
///
/// Each entry in `angles` is `(o_idx, h1, h2, dipole_theta, hh_theta)`.
/// We find the z-bin for the oxygen and append the `(dipole_theta, hh_theta)` pair.
///
/// # Arguments
/// * `atoms`       - Atom container providing positions.
/// * `angles`      - List of (O index, H1, H2, dipole angle, HH angle).
/// * `bins_normal` - z bin edges (length K, right-open).
///
/// # Returns
/// * `Vec<Vec<(f64, f64)>>` where each inner vector holds angle pairs for that z-bin.
pub fn generate_angle_pair_bins_along_normal(
    atoms: &Atoms,
    angles: &Vec<(i32, i32, i32, f64, f64)>, // (o_idx, h1, h2, dipole, hh)
    bins_normal: &[f64],
) -> Vec<Vec<(f64, f64)>> {
    // Pre-allocate one container per z-bin
    let mut angle_bins = vec![ Vec::new(); bins_normal.len() - 1 ];

    if let Some(positions) = &atoms.positions() {
        for &(o_idx, _h1, _h2, dipole_theta, hh_theta) in angles {
            let z = positions[o_idx as usize][2];
            
            // Find z-bin index for this oxygen
            if let Some(bin_idx) = (0..bins_normal.len() - 1)
                .find(|&i| bins_normal[i] <= z && z < bins_normal[i + 1])
            {   
                
                // Push the (dipole, HH) pair to the matched z-bin
                angle_bins[bin_idx].push((dipole_theta, hh_theta));
            }
        }
    }

    angle_bins
}

/// Build a 2D histogram (dipole x HH) restricted to one or more z-ranges.
///
/// For each angle pair `(dipole_theta, hh_theta)` associated with an oxygen at `z`,
/// we include it only if `z` lies in **any** of the provided `z_range` intervals.
/// If included, we locate the (x,y) bins and increment the 2D count.
///
/// # Arguments
/// * `atoms`       - Atom container providing positions.
/// * `angles`      - List of (O index, H1, H2, dipole angle, HH angle).
/// * `z_range`     - List of inclusive-exclusive z-intervals: `[[z_min, z_max], ...]`.
/// * `bins_dipole` - Bin edges for dipole angle.
/// * `bins_hh`     - Bin edges for HH angle.
///
/// # Returns
/// * `Vec<Vec<usize>>` shaped as [ny][nx], where rows map to HH bins (y)
///   and columns map to dipole bins (x).
pub fn generate_angle_pair_bins_for_normal_range(
    atoms: &Atoms,
    angles: &Vec<(i32, i32, i32, f64, f64)>, // (o_idx, h1, h2, dipole, hh),
    z_range: &Vec<[f64; 2]>,                 // [[z_min, z_max], ...]
    bins_dipole: &[f64],                     // e.g., [0.0, 0.1, 0.2, ...]
    bins_hh: &[f64],                         // e.g., [0.0, 0.1, 0.2, ...]
)-> Vec<Vec<usize>> {

    let n_x = bins_dipole.len().saturating_sub(1);
    let n_y = bins_hh.len().saturating_sub(1);

    // Allocate 2D count grid: [y][x] = [HH][dipole]
    let mut hist2d = vec![ vec![0; n_x]; n_y ];

    if let Some(positions) = &atoms.positions() {
        for &(o_idx, _h1, _h2, dip_theta, hh_theta) in angles {
            let z = positions[o_idx as usize][2];

            // Include if z falls within ANY of the provided z-intervals.
            if !z_range.iter().any(|range| range[0] <= z && z < range[1]) {
                continue;
            }

            // Locate dipole bin
            let x_bin = (0..n_x)
                .find(|&i| bins_dipole[i] <= dip_theta && dip_theta < bins_dipole[i+1]);

            // Locate HH bin
            let y_bin = (0..n_y)
                .find(|&j| bins_hh[j] <= hh_theta && hh_theta < bins_hh[j+1]);

            // Increment if both bins exist
            if let (Some(ix), Some(iy)) = (x_bin, y_bin) {
                hist2d[iy][ix] += 1;
            }
        }
    }

    hist2d
}
