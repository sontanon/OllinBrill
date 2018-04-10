#include "tools.h"

#include "param.h"
#include "arrays.h"

#include "geometry_D2alpha.h"

#undef DEBUG

void calculate_alpha_second_covariant_derivative(void)
{
	int k;

	// First calculate conformal laplacian.
	#pragma omp parallel shared(D2alpha_a, D2alpha_b, D2alpha_h, D2alpha_c, D2alpha_lambda, LaplaAlpha)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			// Calculate second covariant derivatve of the lapse.
			calculate_D2alpha(D2alpha_a + k, D2alpha_b + k, D2alpha_h + k,
				D2alpha_c + k, D2alpha_lambda + k, LaplaAlpha + k,
				r[k], Dr_alpha[k], Dz_alpha[k],
				Drr_alpha[k], Dzz_alpha[k], Drz_alpha[k], DDalpha_r[k],
				g_a[k], g_b[k], g_h[k], g_c[k], g_lambda[k],
				Dr_a[k], Dr_b[k], Dr_h[k], Dr_c[k],
				Dz_a[k], Dz_b[k], Dz_h[k], Dz_c[k],
				c[k], h[k], lambda[k], Dz_lambda[k]);
		}
	}

#ifdef DEBUG
	writeSingleFile(LaplaAlpha, "conf_LaplaAlpha.asc");
	writeSingleFile(D2alpha_a, "conf_D2alpha_a.asc");
	writeSingleFile(D2alpha_b, "conf_D2alpha_b.asc");
	writeSingleFile(D2alpha_h, "conf_D2alpha_h.asc");
	writeSingleFile(D2alpha_c, "conf_D2alpha_c.asc");
	writeSingleFile(D2alpha_lambda, "conf_D2alpha_lambda.asc");
#endif

	// Now add terms that complete physical connection.
	#pragma omp parallel shared(D2alpha_a, D2alpha_b, D2alpha_h, D2alpha_c, D2alpha_lambda, LaplaAlpha)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			calculate_phys_D2alpha(D2alpha_a + k, D2alpha_b + k,
				D2alpha_h + k, D2alpha_c + k, D2alpha_lambda + k, LaplaAlpha + k,
				r[k], a[k], b[k], h[k], c[k], lambda[k],
				g_a[k], g_b[k], g_h[k], g_c[k],
				Dr_alpha[k], Dz_alpha[k], Dr_phi[k], Dz_phi[k], psi4[k]);
		}
	}
}
