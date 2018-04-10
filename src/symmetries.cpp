#include "tools.h"

#include "param.h"
#include "arrays.h"

void symmetries(void)
{
	int i, j, k;

	// Loop over number of ghost zones.
	for (k = 0; k < ghost; k++)
	{
		#pragma omp parallel shared(alpha, phi, a, b, h, c, lambda,\
		K, A_a, A_b, A_h, A_c, A_lambda, Deltar, Deltaz,\
		beta_r, beta_z, dtbeta_r, dtbeta_z)
		{
			// Loop in Z direction for axial symmetry.
			#pragma omp for schedule(guided)
			for (j = ghost - k; j < NzTotal; j++)
			{
				alpha[IDX(ghost - 1 - k, j)] = +alpha[IDX(ghost + k, j)];
				phi[IDX(ghost - 1 - k, j)] = +phi[IDX(ghost + k, j)];

				a[IDX(ghost - 1 - k, j)] = +a[IDX(ghost + k, j)];
				b[IDX(ghost - 1 - k, j)] = +b[IDX(ghost + k, j)];
				h[IDX(ghost - 1 - k, j)] = +h[IDX(ghost + k, j)];
				c[IDX(ghost - 1 - k, j)] = +c[IDX(ghost + k, j)];
				lambda[IDX(ghost - 1 - k, j)] = +lambda[IDX(ghost + k, j)];

				K[IDX(ghost - 1 - k, j)] = +K[IDX(ghost + k, j)];
				A_a[IDX(ghost - 1 - k, j)] = +A_a[IDX(ghost + k, j)];
				A_b[IDX(ghost - 1 - k, j)] = +A_b[IDX(ghost + k, j)];
				A_h[IDX(ghost - 1 - k, j)] = +A_h[IDX(ghost + k, j)];
				A_c[IDX(ghost - 1 - k, j)] = +A_c[IDX(ghost + k, j)];
				A_lambda[IDX(ghost - 1 - k, j)] = +A_lambda[IDX(ghost + k, j)];

				Deltar[IDX(ghost - 1 - k, j)] = -Deltar[IDX(ghost + k, j)];
				Deltaz[IDX(ghost - 1 - k, j)] = +Deltaz[IDX(ghost + k, j)];

#ifdef SHIFT
				beta_r[IDX(ghost - 1 - k, j)] = -beta_r[IDX(ghost + k, j)];
				beta_z[IDX(ghost - 1 - k, j)] = +beta_z[IDX(ghost + k, j)];
				dtbeta_r[IDX(ghost - 1 - k, j)] = -dtbeta_r[IDX(ghost + k, j)];
				dtbeta_z[IDX(ghost - 1 - k, j)] = +dtbeta_z[IDX(ghost + k, j)];
#endif
			}
		}

		#pragma omp parallel shared(alpha, phi, a, b, h, c, lambda,\
		K, A_a, A_b, A_h, A_c, A_lambda, Deltar, Deltaz,\
		beta_r, beta_z, dtbeta_r, dtbeta_z)
		{
			// Loop in R direction.
			#pragma omp for schedule(guided)
			for (i = ghost - k; i < NrTotal; i++)
			{
				alpha[IDX(i, ghost - 1 - k)] = +alpha[IDX(i, ghost + k)];
				phi[IDX(i, ghost - 1 - k)] = +phi[IDX(i, ghost + k)];

				a[IDX(i, ghost - 1 - k)] = +a[IDX(i, ghost + k)];
				b[IDX(i, ghost - 1 - k)] = +b[IDX(i, ghost + k)];
				h[IDX(i, ghost - 1 - k)] = +h[IDX(i, ghost + k)];
				c[IDX(i, ghost - 1 - k)] = -c[IDX(i, ghost + k)];
				lambda[IDX(i, ghost - 1 - k)] = +lambda[IDX(i, ghost + k)];

				K[IDX(i, ghost - 1 - k)] = +K[IDX(i, ghost + k)];
				A_a[IDX(i, ghost - 1 - k)] = +A_a[IDX(i, ghost + k)];
				A_b[IDX(i, ghost - 1 - k)] = +A_b[IDX(i, ghost + k)];
				A_h[IDX(i, ghost - 1 - k)] = +A_h[IDX(i, ghost + k)];
				A_c[IDX(i, ghost - 1 - k)] = -A_c[IDX(i, ghost + k)];
				A_lambda[IDX(i, ghost - 1 - k)] = +A_lambda[IDX(i, ghost + k)];

				Deltar[IDX(i, ghost - 1 - k)] = +Deltar[IDX(i, ghost + k)];
				Deltaz[IDX(i, ghost - 1 - k)] = -Deltaz[IDX(i, ghost + k)];

#ifdef SHIFT
				beta_r[IDX(i, ghost - 1 - k)] = +beta_r[IDX(i, ghost + k)];
				beta_z[IDX(i, ghost - 1 - k)] = -beta_z[IDX(i, ghost + k)];
				dtbeta_r[IDX(i, ghost - 1 - k)] = +dtbeta_r[IDX(i, ghost + k)];
				dtbeta_z[IDX(i, ghost - 1 - k)] = -dtbeta_z[IDX(i, ghost + k)];
#endif
			}
		}

		// Corner.
		alpha[IDX(ghost - 1 - k, ghost - 1 - k)] = +alpha[IDX(ghost + k, ghost + k)];
		phi[IDX(ghost - 1 - k, ghost - 1 - k)] = +phi[IDX(ghost + k, ghost + k)];

		a[IDX(ghost - 1 - k, ghost - 1 - k)] = +a[IDX(ghost + k, ghost + k)];
		b[IDX(ghost - 1 - k, ghost - 1 - k)] = +b[IDX(ghost + k, ghost + k)];
		h[IDX(ghost - 1 - k, ghost - 1 - k)] = +h[IDX(ghost + k, ghost + k)];
		c[IDX(ghost - 1 - k, ghost - 1 - k)] = -c[IDX(ghost + k, ghost + k)];
		lambda[IDX(ghost - 1 - k, ghost - 1 - k)] = +lambda[IDX(ghost + k, ghost + k)];

		K[IDX(ghost - 1 - k, ghost - 1 - k)] = +K[IDX(ghost + k, ghost + k)];
		A_a[IDX(ghost - 1 - k, ghost - 1 - k)] = +A_a[IDX(ghost + k, ghost + k)];
		A_b[IDX(ghost - 1 - k, ghost - 1 - k)] = +A_b[IDX(ghost + k, ghost + k)];
		A_h[IDX(ghost - 1 - k, ghost - 1 - k)] = +A_h[IDX(ghost + k, ghost + k)];
		A_c[IDX(ghost - 1 - k, ghost - 1 - k)] = -A_c[IDX(ghost + k, ghost + k)];
		A_lambda[IDX(ghost - 1 - k, ghost - 1 - k)] = +A_lambda[IDX(ghost + k, ghost + k)];

		Deltar[IDX(ghost - 1 - k, ghost - 1 - k)] = -Deltar[IDX(ghost + k, ghost + k)];
		Deltaz[IDX(ghost - 1 - k, ghost - 1 - k)] = -Deltaz[IDX(ghost + k, ghost + k)];

#ifdef SHIFT
		beta_r[IDX(ghost - 1 - k, ghost - 1 - k)] = -beta_r[IDX(ghost + k, ghost + k)];
		beta_z[IDX(ghost - 1 - k, ghost - 1 - k)] = -beta_z[IDX(ghost + k, ghost + k)];
		dtbeta_r[IDX(ghost - 1 - k, ghost - 1 - k)] = -dtbeta_r[IDX(ghost + k, ghost + k)];
		dtbeta_z[IDX(ghost - 1 - k, ghost - 1 - k)] = -dtbeta_z[IDX(ghost + k, ghost + k)];
#endif
	}

	return;
}