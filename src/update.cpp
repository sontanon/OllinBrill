#include "tools.h"

#include "param.h"
#include "arrays.h"

void update(const double dtw)
{
	int k;

	#pragma omp parallel shared(alpha, phi, a, b, h, c, lambda,\
	K, A_a, A_b, A_h, A_c, A_lambda, Deltar, Deltaz)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			alpha[k] = alpha_p[k] + dtw * salpha[k];
			phi[k] = phi_p[k] + dtw * sphi[k];

			a[k] = a_p[k] + dtw * sa[k];
			b[k] = b_p[k] + dtw * sb[k];
			h[k] = h_p[k] + dtw * sh[k];
			c[k] = c_p[k] + dtw * sc[k];
			lambda[k] = lambda_p[k] + dtw * slambda[k];
			K[k] = K_p[k] + dtw * sK[k];
			A_a[k] = A_a_p[k] + dtw * sA_a[k];
			A_b[k] = A_b_p[k] + dtw * sA_b[k];
			A_h[k] = A_h_p[k] + dtw * sA_h[k];
			A_c[k] = A_c_p[k] + dtw * sA_c[k];
			A_lambda[k] = A_lambda_p[k] + dtw * sA_lambda[k];
			Deltar[k] = Deltar_p[k] + dtw * sDeltar[k];
			Deltaz[k] = Deltaz_p[k] + dtw * sDeltaz[k];

#ifdef SHIFT
			// Shift.
			beta_r[k] = beta_r_p[k] + dtw * sbeta_r[k];
			beta_z[k] = beta_z_p[k] + dtw * sbeta_z[k];
			dtbeta_r[k] = dtbeta_r_p[k] + dtw * sdtbeta_r[k];
			dtbeta_z[k] = dtbeta_z_p[k] + dtw * sdtbeta_z[k];
#endif
		}
	}
}