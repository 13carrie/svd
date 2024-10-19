// implement SVD here lol
extern crate nalgebra as na;
use na::{DMatrix, SVD};

pub fn compute_svd() {
    // Example matrix
    let m = DMatrix::from_row_slice(3, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);

    // Compute SVD
    let svd = SVD::new(m.clone(), true, true);

    // Check if SVD computation was successful
    if let (Some(u), Some(v_t)) = (svd.u, svd.v_t) {
        let sigma = svd.singular_values;

        // Print results
        println!("U: \n{}", u);
        println!("Σ: \n{}", sigma);
        println!("V^T: \n{}", v_t);

        println!("SVD computation completed.");
    } else {
        println!("Failed to compute SVD.");
    }
}

fn validate_svd() {
    // validate that svd decomp is correct
}