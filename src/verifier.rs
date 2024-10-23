// implement proof verification logic here
mod circuit;
use circuit::CircuitInput;
use halo2_base::utils::{BigPrimeField, ScalarField};
use axiom_eth::rlp::{
    builder::{FnSynthesize, RlcThreadBuilder, RlpCircuitBuilder},
    rlc::RlcChip,
    *,
};


// polynomial commitment
// check that given commitment is valid for polynomial
// 3 args: pk, c, P
fn verify_commitment() {

}

// polynomial commitment
// check correctness of polynomial evaluation using provided witness
fn verify_open() {

}

// verify proof given proof, verification key and public inputs
// inputs: generated proof, public inputs, verification key
// freivald's algo
// should be really similar to two_step_svd_verif
fn verify_proof<F: ScalarField>(
    mut builder: RlcThreadBuilder<F>,
    input: CircuitInput,
) -> RlpCircuitBuilder<F, impl FnSynthesize<F>> {
    // input a lot more of two_step_svd_verif
}