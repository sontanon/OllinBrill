#include "tools.h"

#include "param.h"
#include "arrays.h"

void accumulate(const int l, const int niter, const double w)
{
	int k;
	if (l == 1)
	{
		// alpha_a = w * salpha.
		#pragma omp parallel shared(phi_a, a_a, b_a, h_a, c_a, lambda_a,\
		A_a_a, A_b_a, A_h_a, A_c_a, A_lambda_a, Deltar_a, Deltaz_a,\
		beta_r_a, beta_z_a, dtbeta_r_a, dtbeta_z_a,\
		K_a, alpha_a)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				phi_a[k] = w * sphi[k];
				a_a[k] = w * sa[k];
				b_a[k] = w * sb[k];
				h_a[k] = w * sh[k];
				c_a[k] = w * sc[k];
				lambda_a[k] = w * slambda[k];
				A_a_a[k] = w * sA_a[k];
				A_b_a[k] = w * sA_b[k];
				A_h_a[k] = w * sA_h[k];
				A_c_a[k] = w * sA_c[k];
				A_lambda_a[k] = w * sA_lambda[k];
				Deltar_a[k] = w * sDeltar[k];
				Deltaz_a[k] = w * sDeltaz[k];

#ifdef SHIFT
				beta_r_a[k] = w * sbeta_r[k];
				beta_z_a[k] = w * sbeta_z[k];
				dtbeta_r_a[k] = w * sdtbeta_r[k];
				dtbeta_z_a[k] = w * sdtbeta_z[k];
#endif

				// Trivial for maximal slicing.
				alpha_a[k] = w * salpha[k];
				K_a[k] = w * sK[k];
			}
		}
	}
	else if (l < niter)
	{
		// alpha_a = alpha_a + w * salpha.
		#pragma omp parallel shared(phi_a, a_a, b_a, h_a, c_a, lambda_a,\
		A_a_a, A_b_a, A_h_a, A_c_a, A_lambda_a, Deltar_a, Deltaz_a,\
		beta_r_a, beta_z_a, dtbeta_r_a, dtbeta_z_a,\
		alpha_a, K_a)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				phi_a[k] += w * sphi[k];
				a_a[k] += w * sa[k];
				b_a[k] += w * sb[k];
				h_a[k] += w * sh[k];
				c_a[k] += w * sc[k];
				lambda_a[k] += w * slambda[k];
				A_a_a[k] += w * sA_a[k];
				A_b_a[k] += w * sA_b[k];
				A_h_a[k] += w * sA_h[k];
				A_c_a[k] += w * sA_c[k];
				A_lambda_a[k] += w * sA_lambda[k];
				Deltar_a[k] += w * sDeltar[k];
				Deltaz_a[k] += w * sDeltaz[k];

#ifdef SHIFT
				beta_r_a[k] += w * sbeta_r[k];
				beta_z_a[k] += w * sbeta_z[k];
				dtbeta_r_a[k] += w * sdtbeta_r[k];
				dtbeta_z_a[k] += w * sdtbeta_z[k];
#endif

				// Trivial for maximal slicing.
				alpha_a[k] += w * salpha[k];
				K_a[k] += w * sK[k];
			}
		}
	}
	else
	{
		// salpha = alpha_a + w * salpha.
		#pragma omp parallel shared(sphi, sa, sb, sh, sc, slambda,\
		sA_a, sA_b, sA_h, sA_c, sA_lambda, sDeltar, sDeltaz,\
		sbeta_r, sbeta_z, sdtbeta_r, sdtbeta_z,\
		salpha, sK)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				sphi[k] = phi_a[k] + w * sphi[k];
				sa[k] = a_a[k] + w * sa[k];
				sb[k] = b_a[k] + w * sb[k];
				sh[k] = h_a[k] + w * sh[k];
				sc[k] = c_a[k] + w * sc[k];
				slambda[k] = lambda_a[k] + w * slambda[k];
				sA_a[k] = A_a_a[k] + w * sA_a[k];
				sA_b[k] = A_b_a[k] + w * sA_b[k];
				sA_h[k] = A_h_a[k] + w * sA_h[k];
				sA_c[k] = A_c_a[k] + w * sA_c[k];
				sA_lambda[k] = A_lambda_a[k] + w * sA_lambda[k];
				sDeltar[k] = Deltar_a[k] + w * sDeltar[k];
				sDeltaz[k] = Deltaz_a[k] + w * sDeltaz[k];

#ifdef SHIFT
				sbeta_r[k] = beta_r_a[k] + w * sbeta_r[k];
				sbeta_z[k] = beta_z_a[k] + w * sbeta_z[k];
				sdtbeta_r[k] = dtbeta_r_a[k] + w * sdtbeta_r[k];
				sdtbeta_z[k] = dtbeta_z_a[k] + w * sdtbeta_z[k];
#endif

				// Trivial for maximal slicing.
				salpha[k] = alpha_a[k] + w * salpha[k];
				sK[k] = K_a[k] + w * sK[k];
			}
		}
	}
}
