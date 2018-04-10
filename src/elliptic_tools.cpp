#include "tools.h"

// Reduce array to elliptic solver-sized array.
// Reduction is done while ghost is equal to its original value,
// i.e. ghost = 2 for second order and ghost = 3 for fourth order.
void ghost_reduce(const double *u, double *g_u, const int NrInterior, const int NzInterior, const int ghost)
{
	int i, j;
	int k = ghost - 1;

#pragma omp parallel shared(g_u) private(j)
	{
#pragma omp for schedule(guided)
		for (i = 0; i < NrInterior + 2; i++)
		{
			for (j = 0; j < NzInterior + 2; j++)
			{
				g_u[i * (NzInterior + 2) + j] = u[IDX(k + i, k + j)];
			}
		}
	}
}

// Fill and array using an elliptic solver-sized array using
// symmetry conditions.
// Fill-in is done after ghost has been reset to its original
// value.
void ghost_fill(const double *g_u, double *u, const int r_sym, const int z_sym, const int NrInterior, const int NzInterior, const int ghost)
{
	int i, j;
	int k = ghost - 1;

	// Fill points that coincide with reduced array.
#pragma omp parallel shared(u) private(j)
	{
#pragma omp for schedule(guided)
		for (i = 0; i < NrInterior + 2; i++)
		{
			for (j = 0; j < NzInterior + 2; j++)
			{
				u[IDX(k + i, k + j)] = g_u[i * (NzInterior + 2) + j];
			}
		}
	}

	// Now fill ghost zones:
	for (k = ghost - 2; k >= 0; k--)
	{
		// Correct R boundaries.
#pragma omp parallel shared(u) 
		{
#pragma omp for schedule(guided)
			for (j = k + 1; j < NzInterior + ghost + 1; j++)
			{
				// Symmetry.
				u[IDX(k, j)] = (double)r_sym * u[IDX(2 * ghost - 1 - k, j)];
			}
		}

		// Correct Z boundaries.
#pragma omp parallel shared(u)
		{
#pragma omp for schedule(guided)
			for (i = k + 1; i < NrInterior + ghost + 1; i++)
			{
				// Symmetry.
				u[IDX(i, k)] = (double)z_sym * u[IDX(i, 2 * ghost - 1 - k)];
			}
		}

		// Correct corner using diagonal symmetry.
		u[IDX(k, k)] = (double)(r_sym * z_sym) * u[IDX(2 * ghost - 1 - k, 2 * ghost - 1 - k)];
	}
}
