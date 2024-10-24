// first, define circuit structure
// then implement/configure circuit for SVD proof constraints and witness generation
use halo2_proofs::{
    arithmetic::FieldExt,
    circuit::{Chip, Layouter, SimpleFloorPlanner},
    plonk::{Circuit, ConstraintSystem, Error},
    poly::{
        Rotation,
        commitment::ParamsProver,
            kzg::{
                commitment::{KZGCommitmentScheme, ParamsKZG},
                multiopen::{ProverSHPLONK, VerifierSHPLONK},
                strategy::SingleStrategy,
            },},
};
use pairing::bn256::Fr as Fp;
use axiom_eth::rlp::{
    builder::{FnSynthesize, RlcThreadBuilder, RlpCircuitBuilder},
    rlc::RlcChip,
    *,
};

mod linear_algebra;
use linear_algebra::{check_svd_phase0, check_svd_phase1};




// zkfixedpointchip for encoding all matrix mult. operations
    // forked from library
// zkmatrix for matrix multiplication and transposition
    // implemented in linear_algebra
// zkvector for inner product, norm, distance, vector operations
    // implemented in linear_algebra

struct CircuitInput {
    // private inputs: U, sigma, V
    u: Vec<Vec<Fp>>,
    sigma: Vec<Fp>,
    v: Vec<Vec<Fp>>,
    // public inputs: matrix that is being decomposed
    m: Vec<Vec<Fp>>,

}

impl circuit<u64> for CircuitInput {
    // configure columns and selectors


    fn configure(meta: &mut ConstraintSystem) -> Self {
        // constraints:
        // U and V must be orthogonal
        // Sigma must be diagonal
        // U*Sigma*V^T must = M
        let u = meta.advice_column();
        let sigma = meta.advice_column();
        let v = meta.advice_column();
        let m = meta.instance_column();

        // orthogonality constraints
        meta.create_gate("orthogonality of u", |meta| {
            let u_i = meta.query_advice(u, Rotation::cur());
            let u_j = meta.query_advice(u, Rotation::next());
            vec![u_i * u_j - Expression::Constant(F::zero())] // For distinct rows/columns
        });

        meta.create_gate("orthogonality of v", |meta| {
            let v_i = meta.query_advice(v, Rotation::cur());
            let v_j = meta.query_advice(v, Rotation::next());
            vec![v_i * v_j - Expression::Constant(F::zero())] // For distinct rows/columns
        });

        // diagonality constraint for Sigma
        meta.create_gate("diagonality of sigma", |meta| {
            let sigma_off_diag = meta.query_advice(sigma, Rotation::cur());
            vec![sigma_off_diag - Expression::Constant(F::zero())] // Off-diagonal elements must be zero
        });

        // U * Sigma * V^T = M constraint
        meta.create_gate("original equation constraint", |meta| {
            let u = meta.query_advice(u, Rotation::cur());
            let sigma = meta.query_advice(sigma, Rotation::cur());
            let v = meta.query_advice(v, Rotation::cur());
            let m = meta.query_instance(m_instance, Rotation::cur());
            vec![u * sigma * v - m] // Ensure the product equals the original matrix
        });

        // initialise struct
        CircuitInput {u, sigma, v, m,}
    }
    
    // 
    fn synthesize(
        &self,
        // building the layout of lookup table
        cs: &mut impl Layouter,
    ) -> Result<(), Error> {
        cs.assign_region(
            || "SVD region",
            |mut region| {
                for (i, &value) in self.u.iter().enumerate() {
                    region.assign_advice(
                        || "u",
                        self.u,
                        i,
                        || Value::known(value),
                    )?;
                }
                for (i, &value) in self.sigma.iter().enumerate() {
                    region.assign_advice(
                        || "sigma",
                        self.sigma,
                        i,
                        || Value::known(value),
                    )?;
                }

                for (i, &value) in self.v.iter().enumerate() {
                    region.assign_advice(
                        || "v",
                        self.v,
                        i,
                        || Value::known(value),
                    )?;
                }

                for (i, &value) in self.m.iter().enumerate() {
                    region.assign_advice(
                        || "m",
                        self.m,
                        i,
                        || Value::known(value),
                    )?;
                }

                // check the "ok" usage here -- correct?
                Ok(())},
            )?;

        check_svd_phase0(cs, &self.u, &self.sigma, &self.v, &self.m)?;
        check_svd_phase1(cs, &self.u, &self.sigma, &self.v, &self.m)?;

        Ok(())
    }
}