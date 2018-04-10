void pardiso_wrapper(const csr_matrix A,	// Matrix system to solve: Au = f.
	double *u,				// Solution array.
	double *f,				// RHS array.
	double *r,				// Residual, r = f - Au, array.
	const double tol,			// Tolerance convergence.
	double *norm,				// Pointer to final norm.
	MKL_INT *convergence,			// Pointer to convergence flag.
	const MKL_INT infnorm,			// Select infnorm or twonorm.
	const MKL_INT perm_use,			// Calculate or use permutation vector:
						// 0: Do not calculate or use permutation.
						// 1: Use specified permutation vector.
						// 2: Calculate permutation vector onto perm.
	const MKL_INT precond_use,		// Use previously computed LU with CGS iteration.
						// 0: Do not use CGS preconditioner.
						// L: Stopping criterion of Krylov-Subspace iteration 10**(-L).
	const MKL_INT low_rank_use);		// Use Low Rank Update.
