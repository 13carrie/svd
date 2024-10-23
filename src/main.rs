mod svd;


fn main() {

    // first, compute svd
    let (m, u, sigma, v_t) = svd::compute_svd();

    // then, setup any necessary env

    
    // then initialise svd proof

    // run svd proof
}

// step 1 of polynomial commitment
// initialises setup params to be used by prover and verifier
/*
fn set_up() -> (i32, i32) {
    let params = 1;
    let verifier_params = 1;
    (params, verifier_params) 
}


fn run_svd_proof() {

}  */