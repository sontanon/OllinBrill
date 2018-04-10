#include "tools.h"
#include "param.h"

#include "flat_laplacian_csr_gen.h"
#include "pardiso_wrapper.h"
#include "elliptic_tools.h"

// Use infinity norm in solver.
#define INFNORM 0

#undef DEBUG

//  Flat Laplacian, solves the linear equation:
//    __2
//  ( \/     + s(r, z) ) u(r, z) = f(r, z).
//  	flat
//
//  Where the Laplacian is the flat Laplacian in cylindrical
//  coordinates:
//
//   __2       2      2
//   \/     = d   +  d    + (1/r) d  .
//     flat    rr     zz           r
//
//  And is solved to a specified order finite difference, 
//  i.e. either second or fourth order.
//
//  s(r, z) is a linear source.
//  f(r, z) is the RHS.
//
//  It returns the solution u(r, z) and the residual res(r, z).
//
void flat_laplacian(double *u,
	const double *f,
	double *res,
	const double *s,
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

	// The main point of this solver is that it works on a smaller grid
	// than that used on the rest of the program.
	// For a second order approximation, we use a grid of NrInterior *
	// NzInterior interior points plus a boundary of one point all arround
	// it, thus a grid of (NrInterior + 2) * (NzInterior + 2).
	// For fourth order, we need an additional boundary point for the 
	// finite difference approximations, thus the grid will be now of
	// (NrInterior + 4) * (NzInterior + 4) points.
	// 
	// The rest of the program runs on grids of (NrInterior + 3) * 
	// (NzInterior + 3) for second order and (NrInterior + 5) *
	// (NzInterior + 5) for fourth order, therefore a reduction is
	// necessary, eliminating the lower-left sides of the grid which
	// are later filled trivially using symmetry conditions.
	//
	// The current value of ghost zones is stored in temporary variable.
	// Ghost zones will later be reset to this original value.
	int temp_ghost = ghost;

	// Size of reduced arrays.
	size_t g_size = (NrInterior + 2) * (NzInterior + 2) * sizeof(double);

	// Allocate reduced arrays.
	double *g_u = (double *)malloc(g_size);
	double *g_f = (double *)malloc(g_size);
	double *g_s = (double *)malloc(g_size);
	double *g_res = (double *)malloc(g_size);

	// Reduce arrays.
	ghost_reduce(u, g_u, NrInterior, NzInterior, ghost);
	ghost_reduce(f, g_f, NrInterior, NzInterior, ghost);
	ghost_reduce(s, g_s, NrInterior, NzInterior, ghost);
	ghost_reduce(res, g_res, NrInterior, NzInterior, ghost);

	// Set new ghost.
	ghost = 1;

	// Set temporary total number of points.
	NrTotal = NrInterior + 2;
	NzTotal = NzInterior + 2;

	// Allocate and generate CSR matrix.
	csr_matrix A;
	int DIM0 = NrTotal * NzTotal;
	int nnz0 = nnz_flat_laplacian(NrInterior, NzInterior, norder, robin);
	csr_allocate(&A, DIM0, DIM0, nnz0);

	// Fill CSR matrix.
	csr_gen_flat_laplacian(A, NrInterior, NzInterior, norder, dr, dz, g_s, g_f, uInf, robin, r_sym, z_sym);

	// Elliptic solver return variables.
	double norm = 0.0;
	int convergence = 0;
	double tol = pow(dr, norder);

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
	NzTotal = NzInterior + ghost + 1;
	NrTotal = NrInterior + ghost + 1;

	// Transfer solution and residual to original arrays.
	ghost_fill(g_u, u, r_sym, z_sym, NrInterior, NzInterior, ghost);
	ghost_fill(g_res, res, r_sym, z_sym, NrInterior, NzInterior, ghost);

	// Clear memory.
	free(g_u);
	free(g_res);
	free(g_f);
	free(g_s);

	// Clear CSR matrix.
	csr_deallocate(&A);
}
