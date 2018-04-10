#include "tools.h"

#include "param.h"
#include "arrays.h"

void array_init(void)
{
	int k;

#pragma omp parallel shared(auxarray, auxarray2, auxarray3, auxarray4, r, z, rr, r2,\
	alpha, Dr_alpha, Dz_alpha, Drr_alpha, Dzz_alpha, Drz_alpha, falpha, a, b, h, c, lambda,\
	phys_a, phys_b, phys_h, phys_c, phys_lambda, g_rr, g_rth, g_thth, g_phph,\
	ig_rr, ig_thth, ig_phph, ig_rth, c_rr, c_th, c_ph, courant_rr, courant_th, courant_ph,\
	Dr_a, Dr_b, Dr_h, Dr_c, Dr_lambda, Dz_a, Dz_b, Dz_h, Dz_c, Dz_lambda,\
	Drr_a, Drr_b, Drr_h, Drr_c, Drr_lambda,\
	Dzz_a, Dzz_b, Dzz_h, Dzz_c, Dzz_lambda,\
	Drz_a, Drz_b, Drz_h, Drz_c, Drz_lambda,\
	g_a, g_b, g_h, g_c, g_lambda,\
	hdet, Dr_hdet, Dz_hdet,\
	Dr_g_a, Dr_g_b, Dr_g_h, Dz_g_a, Dz_g_b, Dz_g_h, Dz_g_c, A_a, A_b, A_h, A_c,\
	A_lambda, Dr_A_a, Dr_A_b, Dr_A_h, Dr_A_c, Dr_A_lambda,\
	Dz_A_a, Dz_A_b, Dz_A_h, Dz_A_c, Dz_A_lambda,\
	Drz_A_a, Drz_A_b, Drz_A_h, Drz_A_c, Drz_K, \
	Drrz_a, Drzz_a, Drrz_b, Drzz_b, Drrz_h, Drzz_h, Drrz_c, Drzz_c, Drrz_phi, Drzz_phi,\
	phi, psi,\
	psi4, Dr_phi, Dz_phi, Drr_phi, Dzz_phi, Drz_phi, K, Dr_K, Dz_K, Deltar, Deltaz,\
	Dr_Deltar, Dr_Deltaz, Dz_Deltar, Dz_Deltaz, R_a, R_b, R_h, R_c, R_lambda, RSCAL,\
	DDeltar_r, DDalpha_r, DDphi_r,\
	D2alpha_a, D2alpha_b, D2alpha_h, D2alpha_c, D2alpha_lambda, LaplaAlpha,\
	A2_a, A2_b, A2_h, A2_c, A2_lambda, A2,\
	upA_a, upA_b, upA_h, upA_c, divA_r, divA_z, ADelta_r, ADelta_z, ham, mom_r, mom_z,\
	Clambda, CA_lambda, CDeltar, CDeltaz,\
	alpha_p, phi_p, a_p, b_p, h_p, c_p, lambda_p, K_p, A_a_p, A_b_p, A_h_p, A_c_p, A_lambda_p, Deltar_p, Deltaz_p,\
	alpha_a, phi_a, a_a, b_a, h_a, c_a, lambda_a, K_a, A_a_a, A_b_a, A_h_a, A_c_a, A_lambda_a, Deltar_a, Deltaz_a,\
	salpha, sphi, sa, sb, sh, sc, slambda, sK, sA_a, sA_b, sA_h, sA_c, sA_lambda, sDeltar, sDeltaz,\
	ell_a, ell_b, ell_c, ell_d, ell_e, ell_s, ell_f,\
	beta_r, beta_z, dtbeta_r, dtbeta_z, Dr_beta_r, Dz_beta_r, Drr_beta_r, Dzz_beta_r, Drz_beta_r,\
	Dr_beta_z, Dz_beta_z, Drr_beta_z, Dzz_beta_z, Drz_beta_z, DD_beta_rr, Dr_dtbeta_r, Dz_dtbeta_r,\
	Dr_dtbeta_z, Dz_dtbeta_z, div_beta, Dr_div_beta, Dz_div_beta, sbeta_r, sbeta_z,\
	sdtbeta_r, sdtbeta_z, DAr_alpha, DAz_alpha, DAr_phi, DAz_phi,\
	DAr_a, DAz_a, DAr_b, DAz_b, DAr_h, DAz_h, DAr_c, DAz_c, DAr_lambda, DAz_lambda,\
	DAr_A_a, DAz_A_a, DAr_A_b, DAz_A_b, DAr_A_h, DAz_A_h, DAr_A_c, DAz_A_c, DAr_A_lambda, DAz_A_lambda,\
	DAr_K, DAz_K, DAr_beta_r, DAz_beta_r, DAr_beta_z, DAz_beta_z, DAr_dtbeta_r, DAz_dtbeta_r,\
	DAr_dtbeta_z, DAz_dtbeta_z, DAr_Deltar, DAz_Deltar, DAr_Deltaz, DAz_Deltaz,\
	beta_r_p, beta_z_p, dtbeta_r_p, dtbeta_z_p, beta_r_a, beta_z_a, dtbeta_r_a, dtbeta_z_a)
	{
#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			// AUXILIARY ARRAYS.
			auxarray[k] = 0.;
			auxarray2[k] = 0.;
			auxarray3[k] = 0.;
			auxarray4[k] = 0.;

			// COORDINATE GRIDS.
			r[k] = 0.;
			z[k] = 0.;
			rr[k] = 0.;
			r2[k] = 0.;

			// LASPE.
			alpha[k] = 0.;
			Dr_alpha[k] = 0.;
			Dz_alpha[k] = 0.;
			Drr_alpha[k] = 0.;
			Dzz_alpha[k] = 0.;
			Drz_alpha[k] = 0.;

			// BONA-MASSO SOURCE.
			falpha[k] = 0.;

			// CONFORMAL METRIC.
			a[k] = 0.;
			b[k] = 0.;
			h[k] = 0.;
			c[k] = 0.;
			lambda[k] = 0.;

			Dr_a[k] = 0.;
			Dr_b[k] = 0.;
			Dr_h[k] = 0.;
			Dr_c[k] = 0.;
			Dr_lambda[k] = 0.;

			Dz_a[k] = 0.;
			Dz_b[k] = 0.;
			Dz_h[k] = 0.;
			Dz_c[k] = 0.;
			Dz_lambda[k] = 0.;

			Drr_a[k] = 0.;
			Drr_b[k] = 0.;
			Drr_h[k] = 0.;
			Drr_c[k] = 0.;
			Drr_lambda[k] = 0.;

			Dzz_a[k] = 0.;
			Dzz_b[k] = 0.;
			Dzz_h[k] = 0.;
			Dzz_lambda[k] = 0.;

			Drz_a[k] = 0.;
			Drz_b[k] = 0.;
			Drz_h[k] = 0.;
			Drz_c[k] = 0.;
			Drz_lambda[k] = 0.;

			// INVERSE CONORMAL METRIC.
			g_a[k] = 0.;
			g_b[k] = 0.;
			g_h[k] = 0.;
			g_c[k] = 0.;
			g_lambda[k] = 0.;

			Dr_g_a[k] = 0.;
			Dr_g_b[k] = 0.;
			Dr_g_h[k] = 0.;
			Dr_g_c[k] = 0.;

			Dz_g_a[k] = 0.;
			Dz_g_b[k] = 0.;
			Dz_g_h[k] = 0.;
			Dz_g_c[k] = 0.;

			hdet[k] = 0.;
			Dr_hdet[k] = 0.;
			Dz_hdet[k] = 0.;

			// TRACELESS EXTRINSIC CURVATURE.
			A_a[k] = 0.;
			A_b[k] = 0.;
			A_h[k] = 0.;
			A_c[k] = 0.;
			A_lambda[k] = 0.;

			Dr_A_a[k] = 0.;
			Dr_A_b[k] = 0.;
			Dr_A_h[k] = 0.;
			Dr_A_c[k] = 0.;
			Dr_A_lambda[k] = 0.;

			Dz_A_a[k] = 0.;
			Dz_A_b[k] = 0.;
			Dz_A_h[k] = 0.;
			Dz_A_c[k] = 0.;
			Dz_A_lambda[k] = 0.;

			// CONFORMAL FACTOR.
			phi[k] = 0.;
			psi[k] = 0.;
			psi4[k] = 0.;
			Dr_phi[k] = 0.;
			Dz_phi[k] = 0.;
			Drr_phi[k] = 0.;
			Dzz_phi[k] = 0.;
			Drz_phi[k] = 0.;

			// TRACE OF EXTRINSIC CURVATURE.
			K[k] = 0.;
			Dr_K[k] = 0.;
			Dz_K[k] = 0.;

			// DELTAS.
			Deltar[k] = 0.;
			Dr_Deltar[k] = 0.;
			Dz_Deltar[k] = 0.;

			Deltaz[k] = 0.;
			Dr_Deltaz[k] = 0.;
			Dz_Deltaz[k] = 0.;

			// REGULARIZATION.
			DDeltar_r[k] = 0.;
			DDalpha_r[k] = 0.;
			DDphi_r[k] = 0.;

			// RICCI TENSOR.
			R_a[k] = 0.;
			R_b[k] = 0.;
			R_h[k] = 0.;
			R_c[k] = 0.;
			R_lambda[k] = 0.;
			RSCAL[k] = 0.;

			// SECOND COVARIANT DERIVATIVE OF THE LAPSE
			D2alpha_a[k] = 0.;
			D2alpha_b[k] = 0.;
			D2alpha_h[k] = 0.;
			D2alpha_c[k] = 0.;
			D2alpha_lambda[k] = 0.;
			LaplaAlpha[k] = 0.;

			// QUADRATIC QUANTITIES OF A.
			A2_a[k] = 0.;
			A2_b[k] = 0.;
			A2_h[k] = 0.;
			A2_c[k] = 0.;
			A2_lambda[k] = 0.;
			A2[k] = 0.;

			// CONTRAVARIANT A.
			upA_a[k] = 0.;
			upA_b[k] = 0.;
			upA_h[k] = 0.;
			upA_c[k] = 0.;

			// DIVERGENCE OF A.
			divA_r[k] = 0.;
			divA_z[k] = 0.;

			// TERM IN DELTA SOURCES.
			ADelta_r[k] = 0.;
			ADelta_z[k] = 0.;

			// CONSTRAINTS.
			ham[k] = 0.;
			mom_r[k] = 0.;
			mom_z[k] = 0.;
			Clambda[k] = 0.;
			CA_lambda[k] = 0.;
			CDeltar[k] = 0.;
			CDeltaz[k] = 0.;

			// OLD ARRAYS.
			alpha_p[k] = 0.;
			phi_p[k] = 0.;
			a_p[k] = 0.;
			b_p[k] = 0.;
			h_p[k] = 0.;
			c_p[k] = 0.;
			lambda_p[k] = 0.;
			K_p[k] = 0.;
			A_a_p[k] = 0.;
			A_b_p[k] = 0.;
			A_h_p[k] = 0.;
			A_c_p[k] = 0.;
			A_lambda_p[k] = 0.;
			Deltar_p[k] = 0.;
			Deltaz_p[k] = 0.;

			// ACCUMULATION ARRAYS.
			alpha_a[k] = 0.;
			phi_a[k] = 0.;
			a_a[k] = 0.;
			b_a[k] = 0.;
			h_a[k] = 0.;
			c_a[k] = 0.;
			lambda_a[k] = 0.;
			K_a[k] = 0.;
			A_a_a[k] = 0.;
			A_b_a[k] = 0.;
			A_h_a[k] = 0.;
			A_c_a[k] = 0.;
			A_lambda_a[k] = 0.;
			Deltar_a[k] = 0.;
			Deltaz_a[k] = 0.;

			// SOURCES.
			salpha[k] = 0.;
			sphi[k] = 0.;
			sa[k] = 0.;
			sb[k] = 0.;
			sh[k] = 0.;
			sc[k] = 0.;
			slambda[k] = 0.;
			sK[k] = 0.;
			sA_a[k] = 0.;
			sA_b[k] = 0.;
			sA_h[k] = 0.;
			sA_c[k] = 0.;
			sA_lambda[k] = 0.;
			sDeltar[k] = 0.;
			sDeltaz[k] = 0.;

			// ELLIPTIC COEFFICIENTS.
			ell_a[k] = 0.;
			ell_b[k] = 0.;
			ell_c[k] = 0.;
			ell_d[k] = 0.;
			ell_e[k] = 0.;
			ell_s[k] = 0.;
			ell_f[k] = 0.;

			// PHYSICAL METRIC.
			phys_a[k] = 0.;
			phys_b[k] = 0.;
			phys_h[k] = 0.;
			phys_c[k] = 0.;
			phys_lambda[k] = 0.;

			// MIXED DERIVATIVES FOR BICUBIC INTERPOLATION.
			Drz_A_a[k] = 0.;
			Drz_A_b[k] = 0.;
			Drz_A_h[k] = 0.;
			Drz_A_c[k] = 0.;
			Drz_K[k] = 0.;
			Drrz_a[k] = 0.;
			Drrz_b[k] = 0.;
			Drrz_h[k] = 0.;
			Drrz_c[k] = 0.;
			Drrz_phi[k] = 0.;
			Drzz_a[k] = 0.;
			Drzz_b[k] = 0.;
			Drzz_h[k] = 0.;
			Drzz_c[k] = 0.;
			Drzz_phi[k] = 0.;

			// PHYSICAL SPHERICAL METRIC.
			g_rr[k] = 0.;
			g_rth[k] = 0.;
			g_thth[k] = 0.;
			g_phph[k] = 0.;
			ig_rr[k] = 0.;
			ig_rth[k] = 0.;
			ig_thth[k] = 0.;
			ig_phph[k] = 0.;

			// SPEED OF LIGHT.
			c_rr[k] = 0.;
			c_th[k] = 0.;
			c_ph[k] = 0.;

			// COURANT CONDITION.
			courant_rr[k] = 0.;
			courant_th[k] = 0.;
			courant_ph[k] = 0.;

			// SHIFT.
			beta_r[k] = 0.;
			beta_z[k] = 0.;
			dtbeta_r[k] = 0.;
			dtbeta_z[k] = 0.;
			Dr_beta_r[k] = 0.;
			Dz_beta_r[k] = 0.;
			Drr_beta_r[k] = 0.;
			Dzz_beta_r[k] = 0.;
			Drz_beta_r[k] = 0.;
			Dr_beta_z[k] = 0.;
			Dz_beta_z[k] = 0.;
			Drr_beta_z[k] = 0.;
			Dzz_beta_z[k] = 0.;
			Drz_beta_z[k] = 0.;
			DD_beta_rr[k] = 0.;
			Dr_dtbeta_r[k] = 0.;
			Dz_dtbeta_r[k] = 0.;
			Dr_dtbeta_z[k] = 0.;
			Dz_dtbeta_z[k] = 0.;
			div_beta[k] = 0.;
			Dr_div_beta[k] = 0.;
			Dz_div_beta[k] = 0.;
			sbeta_r[k] = 0.;
			sbeta_z[k] = 0.;
			sdtbeta_r[k] = 0.;
			sdtbeta_z[k] = 0.;
			beta_r_p[k] = 0.;
			beta_z_p[k] = 0.;
			dtbeta_r_p[k] = 0.;
			dtbeta_z_p[k] = 0.;
			beta_r_a[k] = 0.;
			beta_z_a[k] = 0.;
			dtbeta_r_a[k] = 0.;
			dtbeta_z_a[k] = 0.;

			// ADVECTIVE DERIVATIVES.
			DAr_alpha[k] = 0.;
			DAz_alpha[k] = 0.;
			DAr_phi[k] = 0.;
			DAz_phi[k] = 0.;
			DAr_a[k] = 0.;
			DAz_a[k] = 0.;
			DAr_b[k] = 0.;
			DAz_b[k] = 0.;
			DAr_h[k] = 0.;
			DAz_h[k] = 0.;
			DAr_c[k] = 0.;
			DAz_c[k] = 0.;
			DAr_lambda[k] = 0.;
			DAz_lambda[k] = 0.;
			DAr_A_lambda[k] = 0.;
			DAz_A_lambda[k] = 0.;
			DAr_beta_r[k] = 0.;
			DAz_beta_r[k] = 0.;
			DAr_beta_z[k] = 0.;
			DAz_beta_z[k] = 0.;
			DAr_dtbeta_r[k] = 0.;
			DAz_dtbeta_r[k] = 0.;
			DAr_dtbeta_z[k] = 0.;
			DAz_dtbeta_z[k] = 0.;
			DAr_K[k] = 0.;
			DAz_K[k] = 0.;
			DAr_A_a[k] = 0.;
			DAz_A_a[k] = 0.;
			DAr_A_b[k] = 0.;
			DAz_A_b[k] = 0.;
			DAr_A_h[k] = 0.;
			DAz_A_h[k] = 0.;
			DAr_A_c[k] = 0.;
			DAz_A_c[k] = 0.;
			DAr_Deltar[k] = 0.;
			DAz_Deltar[k] = 0.;
			DAr_Deltaz[k] = 0.;
			DAz_Deltaz[k] = 0.;
		}
	}

	return;
}
