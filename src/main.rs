mod svd;


fn main() {

    // first, compute svd
    let (m, u, sigma, v_t) = svd::compute_svd();

    // establish polynomial commitments
    // then initialise svd proof


    // run svd proof
}

// step 1 of polynomial commitment
// initialises setup params to be used by prover and verifier
/*
fn set_up() -> (?,?) {
    let params = ?;
    let verifier_params = ?;
    (params, verifier_params) 
}


fn run_svd_proof() {

}  */