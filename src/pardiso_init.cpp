#define PARDISO_MAIN_FILE
#include "pardiso_param.h"

void pardiso_init(const MKL_INT NrInterior, const MKL_INT NzInterior)
{
	// Real unsymmetric matrix.
	mtype = 11;
	// One RHS.
	nrhs = 1;
	// Matrix dimension.
	n = (NrInterior + 2) * (NzInterior + 2);

	// Setup PARDISO control parameters. First set all to zero.
	for (int i = 0; i < 64; i++)
		iparm[i] = 0;

	// Default fine-tune parameters. Used for perm_use = 0.
	iparm[1 - 1] = 1;	// Do not use default parameters.
	iparm[2 - 1] = 3;	// Parallel fill-in reordering from METIS.
	iparm[4 - 1] = 0;	// No iterative-direct algorithm.
	iparm[5 - 1] = 0;	// No user fill-in reducing permutation.
	iparm[6 - 1] = 0;	// Do not write solution into RHS.
	iparm[7 - 1] = 0;	// Not in use.
	iparm[8 - 1] = 10;	// Max numbers of iterative refinement steps.
	iparm[9 - 1] = 0;	// Not in use.
	iparm[10 - 1] = 13;	// Perturb the pivot elements with 1E-13.
	iparm[11 - 1] = 1;	// Use nonsymmetric permutation and scaling MPS.
	iparm[12 - 1] = 0;	// No conjugate transposed/transpose solve.
	iparm[13 - 1] = 1;      // Maximum weighted matching algorithm is switched-on (default for non-symmetric).
	iparm[14 - 1] = 0;	// Output: Number of perturbed pivots.
	iparm[15 - 1] = 0;	// Not in use.
	iparm[16 - 1] = 0;	// Not in use.
	iparm[17 - 1] = 0;	// Not in use.
	iparm[18 - 1] = 0;	// No Output: Number of nonzeros in the factor LU.
	iparm[19 - 1] = 0;	// No Output: Mflops for LU factorization.
	iparm[20 - 1] = 0;      // Output: Numbers of CG Iterations.
	iparm[24 - 1] = 10;	// Parallel Numerical Factorization.
	iparm[25 - 1] = 1;	// Parallel Forward/Backward Solve.
	iparm[27 - 1] = 0;	// No matrix check.
	maxfct = 1;		// Maximum number of numerical factorizations.
	mnum = 1;		// Which factorization to use.
	msglvl = MESSAGE_LEVEL;	// Print statistical information in file.
	error = 0;		// Initialize error flag.

	// Initialize internal solver memory pointer.
	for (int i = 0; i < 64; i++)
		pt[i] = 0;

	// Allocate permutation vector.
	perm = (int *)malloc(n * sizeof(int));

#ifdef VERBOSE
	printf("PARDISO: Setup solver memory and parameters.\n");
#endif

	// Setup matrix-vector multiplication type.
	uplo = "non-transposed";


	return;
}
