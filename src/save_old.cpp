#include "tools.h"

#include "param.h"
#include "arrays.h"

void save_old(void)
{
	int k;

	#pragma omp parallel shared(alpha_p, phi_p, a_p, b_p, h_p, c_p, lambda_p, \
	K_p, A_a_p, A_b_p, A_h_p, A_c_p, A_lambda_p, Deltar_p, Deltaz_p,\
	beta_r_p, beta_z_p, dtbeta_r_p, dtbeta_z_p)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			alpha_p[k] = alpha[k];
			phi_p[k] = phi[k];
			a_p[k] = a[k];
			b_p[k] = b[k];
			h_p[k] = h[k];
			c_p[k] = c[k];
			lambda_p[k] = lambda[k];
			K_p[k] = K[k];
			A_a_p[k] = A_a[k];
			A_b_p[k] = A_b[k];
			A_h_p[k] = A_h[k];
			A_c_p[k] = A_c[k];
			A_lambda_p[k] = A_lambda[k];
			Deltar_p[k] = Deltar[k];
			Deltaz_p[k] = Deltaz[k];
#ifdef SHIFT
			beta_r_p[k] = beta_r[k];
			beta_z_p[k] = beta_z[k];
			dtbeta_r_p[k] = dtbeta_r[k];
			dtbeta_z_p[k] = dtbeta_z[k];
#endif
		}
	}
}
