#include "tools.h"
#include "param.h"

#include "general_elliptic_csr_gen.h"
#include "pardiso_wrapper.h"
#include "elliptic_tools.h"

// Use infinity norm in solver.
#define INFNORM 0

void general_elliptic(double *u,
	double *res,
	const double *ell_a,
	const double *ell_b,
	const double *ell_c,
	const double *ell_d,
	const double *ell_e,
	const double *ell_s,
	const double *ell_f,
	const double uInf,
	const int robin,
	const int r_sym,
	const int z_sym,
	const int perm_use,
	const int precond_use,
	const int low_rank_use)
{
	int norder = 0;

	if (strcmp(order, "two") == 0)
	{
		norder = 2;
	}
	else if (strcmp(order, "four") == 0)
	{
		norder = 4;
	}

	// The current value of ghost zones is stored in temporary variable.
	// Ghost zones will later be reset to this original value.
	int temp_ghost = ghost;

	// Size of reduced arrays.
	size_t g_size;

	// Set array size.
	g_size = (NrInterior + 2) * (NzInterior + 2) * sizeof(double);

	// Allocate reduced arrays.
	double *g_u = (double *)malloc(g_size);
	double *g_f = (double *)malloc(g_size);
	double *g_a = (double *)malloc(g_size);
	double *g_b = (double *)malloc(g_size);
	double *g_c = (double *)malloc(g_size);
	double *g_d = (double *)malloc(g_size);
	double *g_e = (double *)malloc(g_size);
	double *g_s = (double *)malloc(g_size);
	double *g_res = (double *)malloc(g_size);

	// Reduce arrays.
	ghost_reduce(u, g_u, NrInterior, NzInterior, ghost);
	ghost_reduce(res, g_res, NrInterior, NzInterior, ghost);
	ghost_reduce(ell_a, g_a, NrInterior, NzInterior, ghost);
	ghost_reduce(ell_b, g_b, NrInterior, NzInterior, ghost);
	ghost_reduce(ell_c, g_c, NrInterior, NzInterior, ghost);
	ghost_reduce(ell_d, g_d, NrInterior, NzInterior, ghost);
	ghost_reduce(ell_e, g_e, NrInterior, NzInterior, ghost);
	ghost_reduce(ell_s, g_s, NrInterior, NzInterior, ghost);
	ghost_reduce(ell_f, g_f, NrInterior, NzInterior, ghost);

	// Set new ghost.
	ghost = 1;

	// Set temporary total number of points.
	NrTotal = NrInterior + 2 * ghost;
	NzTotal = NzInterior + 2 * ghost;

	// Allocate and generate CSR matrix.
	csr_matrix A;
	int DIM0 = NrTotal * NzTotal;
	int nnz0 = nnz_general_elliptic(NrInterior, NzInterior, norder, robin);
	csr_allocate(&A, DIM0, DIM0, nnz0);

	// Fill CSR matrix.
	csr_gen_general_elliptic(A, NrInterior, NzInterior, norder, dr, dz, g_a, g_b, g_c, g_d, g_e, g_s, g_f, uInf, robin, r_sym, z_sym);

	//printf("GENERATED CSR MATRIX.\n");

	// Elliptic solver return variables.
	double norm = 0.0;
	int convergence = 0;
	double tol = pow(dr, 4);

	// Call elliptic solver.
	pardiso_wrapper(A, g_u, g_f, g_res, tol, &norm, &convergence, INFNORM, perm_use, precond_use, low_rank_use);

	// Check solver convergence.
	if (convergence == 1)
	{
		printf("***                                ***\n");
		printf("***   Elliptic solver converged    ***\n");
		printf("***                                ***\n");
	}
	else
	{
		printf("***                                ***\n");
		printf("***   Elliptic solver failed!  %d   ***\n", convergence);
		printf("***                                ***\n");
	}
	printf("***   Residual norm = %3.3e    ***\n", norm);
	printf("***                                ***\n");
	printf("**************************************\n");

	// Reset ghost and total number of points.
	ghost = temp_ghost;
	NrTotal = NrInterior + ghost + 1;
	NzTotal = NzInterior + ghost + 1;

	// Transfer solution and residual to original arrays.
	ghost_fill(g_u, u, r_sym, z_sym, NrInterior, NzInterior, ghost);
	ghost_fill(g_res, res, r_sym, z_sym, NrInterior, NzInterior, ghost);

	// Clear memory.
	free(g_u);
	free(g_res);
	free(g_a);
	free(g_b);
	free(g_c);
	free(g_d);
	free(g_e);
	free(g_s);
	free(g_f);

	// Clear CSR matrix.
	csr_deallocate(&A);
}
