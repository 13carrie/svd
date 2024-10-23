// This file implements essential linear algebra functionalities.
// Much of it has been directly copied from halo2-svd

use halo2_base::gates::{GateChip, GateInstructions, RangeChip, RangeInstructions};
use halo2_base::utils::{biguint_to_fe, BigPrimeField};
use halo2_base::{AssignedValue, QuantumCell};
use halo2_base::{
    Context,
    QuantumCell::{Constant, Existing},
};

// allows ZkVector to be cloned
// holds a vector of AssignedValue<F>
// where F = type that implements BigPrimeField and PRECISION_BITS = constant
// representing precision of fixed-point arithmetic
#[derive(Clone)]
pub struct ZkVector<F: BigPrimeField, const PRECISION_BITS: u32> {
    pub v: Vec<AssignedValue<F>>,
}

// implementation of ZkVector struct
// essentially copy-pasted from halo2-svd
impl<F: BigPrimeField, const PRECISION_BITS: u32> ZkVector<F, PRECISION_BITS> {
    // new
    pub fn new(
        ctx: &mut Context<F>,
        fpchip: &FixedPointChip<F, PRECISION_BITS>,
        v: &Vec<f64>,
    ) -> Self {
        let mut zk_v: Vec<AssignedValue<F>> = Vec::new();
        for elem in v {
            let elem = fpchip.quantization(*elem);
            zk_v.push(ctx.load_witness(elem));
        }
        return Self { v: zk_v };
    }

    // size
    pub fn size(&self) -> usize {
        return self.v.len();
    }

    // dequantize
    // converting fixed-point representation back into the real number that it represents
    pub fn dequantize(&self, fpchip: &FixedPointChip<F, PRECISION_BITS>) -> Vec<f64> {
        let mut dq_v: Vec<f64> = Vec::new();
        for elem in &self.v {
            dq_v.push(fpchip.dequantization(*elem.value()));
        }
        return dq_v;
    }

    /// Prints the dequantized version of the vector and returns it;
    ///
    /// Action is not constrained in anyway
    pub fn print(&self, fpchip: &FixedPointChip<F, PRECISION_BITS>) {
        let dq_v = self.dequantize(fpchip);
        println!("[");
        for elem in dq_v {
            println!("{:?}, ", elem);
        }
        println!("]");
    }

    /// With zk constraints calculates the inner product of this vector with vector x
    ///
    /// Outputs the inner product
    ///
    /// Order doesn't matter because we are only dealing with real numbers
    ///
    /// Low level function; uses the fact that FixedPointChip.{add, mul} just call GateChip.{add, mul}
    ///
    /// Leads to about [self.size()] + 90 constraints
    pub fn inner_product(
        &self,
        ctx: &mut Context<F>,
        fpchip: &FixedPointChip<F, PRECISION_BITS>,
        x: &Vec<AssignedValue<F>>,
    ) -> AssignedValue<F> {
        // couldn't figure out how to use inner_product of fpchip because we use x: &Vec and I didn't want to move
        assert!(self.size() == x.len());

        let mut v: Vec<QuantumCell<F>> = Vec::new();
        for elem in &self.v {
            v.push(Existing(*elem));
        }
        let v = v;

        let mut u: Vec<QuantumCell<F>> = Vec::new();
        for elem in x {
            u.push(Existing(*elem));
        }
        let u = u;

        let res_s = fpchip.gate().inner_product(ctx, u, v);

        // #CONSTRAINTS = 90
        // Implementing this way allows us to amortize the cost of calling this expensive rescaling- will also lead to more accuracy
        let (res, _) = fpchip.signed_div_scale(ctx, res_s);
        return res;
    }

    /// With zk constraints calculates the square of the norm of the vector;
    ///
    /// Outputs the square of the norm
    pub fn _norm_square(
        &self,
        ctx: &mut Context<F>,
        fpchip: &FixedPointChip<F, PRECISION_BITS>,
    ) -> AssignedValue<F> {
        return self.inner_product(ctx, fpchip, &self.v);
    }

    /// With zk constraints calculates the norm of the vector;
    ///
    /// Outputs the norm;
    ///
    /// Square root function is expensive and adds a lot error; Avoid using whenever possible
    pub fn norm(
        &self,
        ctx: &mut Context<F>,
        fpchip: &FixedPointChip<F, PRECISION_BITS>,
    ) -> AssignedValue<F> {
        let norm_sq = self._norm_square(ctx, fpchip);
        return fpchip.qsqrt(ctx, norm_sq);
    }

    /// With zk constraints calculates the distance squared of the vector from vector `x`;
    ///
    /// Outputs the distance squared
    pub fn _dist_square(
        &self,
        ctx: &mut Context<F>,
        fpchip: &FixedPointChip<F, PRECISION_BITS>,
        x: &Vec<AssignedValue<F>>,
    ) -> AssignedValue<F> {
        assert_eq!(self.size(), x.len());
        let mut diff: Vec<AssignedValue<F>> = Vec::new();
        for (r, s) in self.v.iter().zip(x.iter()) {
            diff.push(fpchip.qsub(ctx, *r, *s));
        }
        let diff = Self { v: diff };
        return diff._norm_square(ctx, fpchip);
    }

    /// With zk constraints calculates the dist of the vector from vector `x`
    ///
    /// Outputs the norm;
    ///
    /// Square root function adds a lot error and is very expensive; Avoid using whenever possible
    pub fn dist(
        &self,
        ctx: &mut Context<F>,
        fpchip: &FixedPointChip<F, PRECISION_BITS>,
        x: &Vec<AssignedValue<F>>,
    ) -> AssignedValue<F> {
        let dist_sq = self._dist_square(ctx, fpchip, x);
        return fpchip.qsqrt(ctx, dist_sq);
    }

    /// Multiplies this vector by matrix `a` in the zk-circuit and returns the constrained output `a.v`
    ///
    /// Adds about N^2+90*N constraints
    pub fn mul(
        &self,
        ctx: &mut Context<F>,
        fpchip: &FixedPointChip<F, PRECISION_BITS>,
        a: &ZkMatrix<F, PRECISION_BITS>,
    ) -> Self {
        assert_eq!(a.num_col, self.size());
        let mut y: Vec<AssignedValue<F>> = Vec::new();
        // #CONSTRAINTS = N^2+90*N
        for row in &a.matrix {
            y.push(self.inner_product(ctx, fpchip, row));
        }
        return Self { v: y };
    }

    /// Constrains all the entries of the vector to be in [0, 2^max_bits)
    pub fn entries_less_than(
        &self,
        ctx: &mut Context<F>,
        fpchip: &FixedPointChip<F, PRECISION_BITS>,
        max_bits: usize,
    ) {
        for elem in &self.v {
            fpchip.range_gate().range_check(ctx, *elem, max_bits);
        }
    }

    /// Assumes all entries of the vector are in [0, 2^max_bits) (fails silently otherwise)
    ///
    /// Constrains the entries to be in decreasing order
    pub fn entries_in_desc_order(
        &self,
        ctx: &mut Context<F>,
        fpchip: &FixedPointChip<F, PRECISION_BITS>,
        max_bits: usize,
    ) {
        let mut vec_diff: Vec<AssignedValue<F>> = Vec::new();

        for i in 0..(self.v.len() - 1) {
            let ele = fpchip.qsub(ctx, self.v[i], self.v[i + 1]);
            vec_diff.push(ele);
        }

        for elem in &vec_diff {
            fpchip.range_gate().range_check(ctx, *elem, max_bits);
        }
    }
}

#[derive(Clone)]
pub struct ZkMatrix {
    // Fields and methods for ZkMatrix
}

impl ZkMatrix {
    // Implement methods like matrix multiplication, transpose, etc.
}
