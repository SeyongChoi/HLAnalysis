use pyo3::prelude::*;
/// Compute a simple 1D histogram of the given data using specified bin edges.
///
/// # Arguments
/// * `data` - Vector of floating-point values to be binned.
/// * `bins` - Vector of bin edges (must be sorted in ascending order).
///
/// # Returns
/// * `(counts, bins)` -
///   - `counts`: Number of values falling into each bin.
///   - `bins`: The original bin edges (returned for convenience).
#[pyfunction]
pub fn histogram_1d(data: Vec<f64>, bins: Vec<f64>) -> (Vec<usize>, Vec<f64>) {
    
    // Initialize count vector with one entry per bin (number of bins = bins.len() - 1)
    let mut counts = vec![0; bins.len() - 1];

    // Iterate over all data points
    for value in data {
        // Find which bin this value belongs to
        for i in 0..bins.len() - 1 {
            // If the value lies between two consecutive bin edges, increment the count
            if value >= bins[i] && value < bins[i + 1] {
                counts[i] += 1;
                break; // Once placed, move to the next value
            }
        }
    }

    // Return the bin counts and the bin edges
    (counts, bins)
}
