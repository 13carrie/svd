This is an implementation of SVD (singular value decomposition)
in zero-knowledge, utilising the Rust crate 'halo2'.

Objectives: verify the correctness of SVD without revealing
 the matrix or its decomposition. Prove that the prover knows
 some private inputs that make SVD hold.

Relation R: specifies which public/private input combos are valid
Circuit: implementation of relation R
Arithmetisation: language used to express a circuit for a proof system
Advice valeus: intermediate values that can be efficiently
 computed from private and public inputs
Witness: private inputs + advice values, collectively
SVD algorithm: U*sigma*V
U: columns of orthonormal left singular vectors
Sigma: diagonal matrix of singular values
V*: the 'conjugate transpose' of V, an n*n unitary matrix such that
 V* rows are orthonormal right singular vectors

To prove that result of SVD is valid:
Naive method: compute U*sigma*V and verify that it equals M
 (too intensive)
Better method:
M, U and V can ultimately be broken down into matrix
 multiplications. Use an efficient algorithm for verifying
 matrix multiplications, and then use it to verify the M, U and V
 multiplications

U is an element of the set of real orthogonal matrices of size n x n
V is an element of the set of real orthogonal matrices of size m x m
M is an element of the set of real matrices of size m x n
Sigma is an element of the set of real matrices of size m x n

Observation: given V*(V^T) = Id
