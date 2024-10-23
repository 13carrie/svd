// first, define circuit structure
// then implement/configure circuit for SVD proof constraints and witness generation
use halo2_proofs::{
    arithmetic::FieldExt,
    circuit::{Chip, Layouter, SimpleFloorPlanner},
    plonk::{Circuit, ConstraintSystem, Error},
    poly::Rotation,
};
use pairing::bn256::Fr as Fp;



// zkfixedpointchip for encoding all matrix mult. operations
    // forked from library
// zkmatrix for matrix multiplication and transposition
    // implemented in linear_algebra
// zkvector for inner product, norm, distance, vector operations
    // implemented in linear_algebra

struct circuit {
    // private inputs: U, sigma, V
    u: Vec<Vec<Fp>>,
    sigma: Vec<Fp>,
    v: Vec<Vec<Fp>>,
    // public inputs: matrix that is being decomposed
    m: Vec<Vec<Fp>>,

}


impl circuit<u64> for circuit {
    // configure columns and selectors

    fn define_constraints() {
        // constraints:
        // U and V must be orthogonal
        // Sigma must be diagonal
        // U*Sigma*V^T must = M
    }
    
    fn generate_witness() {
    
    }
}