#include "tools.h"
#include "pardiso_param.h"

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
	const MKL_INT low_rank_use)		// Use low-rank update.
{
	// Auxiliary doubles.
	double res, res0;

	// Modify parameters according to permutation use.
	switch (perm_use)
	{
		// Do not use permutation vector and do not calculate it either.
		case 0:
			iparm[5 - 1] = 0;
			iparm[31 - 1] = 0;
			iparm[36 - 1] = 0;
			break;
		// Use specified permutation, iparm(2) is ignored.
		case 1:
			iparm[5 - 1] = 1;
			iparm[31 - 1] = 0;
			iparm[36 - 1] = 0;
			break;
		// Calculate permutation.
		// Permutation vector computed in phase 11 is returned in the perm array.
		case 2:
			iparm[5 - 1] = 2;
			iparm[31 - 1] = 0;
			iparm[36 - 1] = 0;
			break;
	}

	// Modify parameters according to CGS preconditioner.
	if (precond_use)
	{
		/// LU preconditioned with CGS.
		iparm[4 - 1] = 10 * precond_use + 1;
	}

	// Modify paramters according to Low Rank Update.
	if (low_rank_use)
	{
		iparm[39 - 1] = 1;
		iparm[24 - 1] = 10;
		// No permutation restriction.
		iparm[5 - 1] = 0;
		iparm[31 - 1] = 0;
		iparm[36 - 1] = 0;
		// Aditional default values.
		iparm[4 - 1] = 0;
		iparm[6 - 1] = 0;
		iparm[28 - 1] = 0;
		iparm[37 - 1] = 0;
		iparm[56 - 1] = 0;
		iparm[60 - 1] = 0;
	}
	
	// PARDISO calls are different when using Low Rank Update.
	if (!low_rank_use)
	{
		// Reordering and symbolic factorization.
		phase = 11;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, A.a, A.ia, A.ja, perm, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

		if (error != 0) 
		{
			printf("ERROR during symbolic factorization: %d.\n", error);
			exit(1);
		}
		
#ifdef VERBOSE
		printf("PARDISO: Reordering completed.\n");
		printf("PARDISO: Number of nonzeros in factors = %d.\n", iparm[18 - 1]);
		printf("PARDISO: Number of factorization MFLOPS = %d.\n", iparm[19 - 1]);
#endif

		// Numerical factorization.
		phase = 22;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, A.a, A.ia, A.ja, perm, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

		if (error != 0) 
		{
			printf("ERROR during numerical factorization: %d.\n", error);
			exit(2);
		}

#ifdef VERBOSE
		printf("PARDISO: Factorization completed.\n");
#endif

		// Back substitution and iterative refinement.
		phase = 33;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, A.a, A.ia, A.ja, perm, &nrhs, iparm, &msglvl, f, u, &error);

		// Report CGS iterations.
#ifdef VERBOSE
		if (precond_use)
			printf("PARDISO CGS PRECONDITIONER: iparm(20) = %d.\n", iparm[20 - 1]);
#endif

		if (error != 0) 
		{
			printf("ERROR during solution: %d,\n", error);
			exit(3);
		}
	}
	// Use Low Rank Update by calling directly onto the factorization phase.
	else 
	{
		// Numerical factorization.
		phase = 22;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, A.a, A.ia, A.ja, diff, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

		if (error != 0) 
		{
			printf("ERROR during numerical factorization: %d.\n", error);
			exit(2);
		}

#ifdef VERBOSE
		printf("PARDISO: Factorization completed.\n");
#endif

		// Back substitution and iterative refinement.
		phase = 33;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, A.a, A.ia, A.ja, diff, &nrhs, iparm, &msglvl, f, u, &error);

		if (error != 0) 
		{
			printf("ERROR during solution: %d,\n", error);
			exit(3);
		}
	}

	// Compute residual
	mkl_dcsrgemv(uplo, &A.nrows, A.a, A.ia, A.ja, u, r);
	cblas_daxpy(A.nrows, -1.0, f, 1, r, 1);


	// Calculate norms.
	if (infnorm) 
	{
		res = ABS(r[cblas_idamax(A.nrows, r, 1)]);
		res0 = ABS(f[cblas_idamax(A.nrows, f, 1)]);
	}
	else 
	{
		res = cblas_dnrm2(A.nrows, r, 1);
		res0 = cblas_dnrm2(A.nrows, f, 1);
	}
	// Relative residual.
	res0 = res / res0;

#ifdef VERBOSE
	printf("PARDISO: Relative residual = %e.\n", res0);
	printf("PARDISO: Absolute residual = %e.\n", res);
#endif

	// Check residual and output convergence type.
	if (res < tol) 
	{
		if (res < res0) 
		{
			// Absolute and relative convergence.
			*norm = res;
			*convergence = 1;
#ifdef VERBOSE
			printf("PARDISO: Converged both relatively and absolutely.\n");
#endif
		}
		else 
		{
			// Only absolute convergence.
			*norm = res;
			*convergence = 1;
#ifdef VERBOSE
			printf("PARDISO: Converged absolutely but not relatively.\n");
#endif
		}
	}
	else if (res0 < tol) 
	{
		// Only relative convergence.
		*norm = res;
		*convergence = 2;
#ifdef VERBOSE
		printf("PARDISO: Converged relatively but not absolutely.\n");
#endif
	}
	else 
	{
		// No convergence.
		*norm = res;
		*convergence = 0;
#ifdef VERBOSE
		printf("\nPARDISO: WARNING: Failed to converge!\n\n");
#endif
	}

	// Return.
	return;
}
