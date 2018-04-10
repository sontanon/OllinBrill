#include "tools.h"

#include "arrays.h"
#include "param.h"

void allocate_arrays(const bool status)
{
	// Size of arrays.
	size_t size = DIM * sizeof(double);

	// Allocate arrays.
	if (status)
	{
		// AUXILIARY ARRAYS.
		auxarray = (double *)malloc(size);
		auxarray2 = (double *)malloc(size);
		auxarray3 = (double *)malloc(size);
		auxarray4 = (double *)malloc(size);

		// COORDINATE GRIDS.
		r = (double *)malloc(size);
		z = (double *)malloc(size);
		rr = (double *)malloc(size);
		r2 = (double *)malloc(size);

		// LASPE.
		alpha = (double *)malloc(size);
		Dr_alpha = (double *)malloc(size);
		Dz_alpha = (double *)malloc(size);
		Drr_alpha = (double *)malloc(size);
		Dzz_alpha = (double *)malloc(size);
		Drz_alpha = (double *)malloc(size);

		// BONA-MASSO SOURCE.
		falpha = (double *)malloc(size);

		// CONFORMAL METRIC.
		a = (double *)malloc(size);
		b = (double *)malloc(size);
		h = (double *)malloc(size);
		c = (double *)malloc(size);
		lambda = (double *)malloc(size);

		Dr_a = (double *)malloc(size);
		Dr_b = (double *)malloc(size);
		Dr_h = (double *)malloc(size);
		Dr_c = (double *)malloc(size);
		Dr_lambda = (double *)malloc(size);

		Dz_a = (double *)malloc(size);
		Dz_b = (double *)malloc(size);
		Dz_h = (double *)malloc(size);
		Dz_c = (double *)malloc(size);
		Dz_lambda = (double *)malloc(size);

		Drr_a = (double *)malloc(size);
		Drr_b = (double *)malloc(size);
		Drr_h = (double *)malloc(size);
		Drr_c = (double *)malloc(size);
		Drr_lambda = (double *)malloc(size);

		Dzz_a = (double *)malloc(size);
		Dzz_b = (double *)malloc(size);
		Dzz_h = (double *)malloc(size);
		Dzz_c = (double *)malloc(size);
		Dzz_lambda = (double *)malloc(size);

		Drz_a = (double *)malloc(size);
		Drz_b = (double *)malloc(size);
		Drz_h = (double *)malloc(size);
		Drz_c = (double *)malloc(size);
		Drz_lambda = (double *)malloc(size);

		// INVERSE CONORMAL METRIC.
		g_a = (double *)malloc(size);
		g_b = (double *)malloc(size);
		g_h = (double *)malloc(size);
		g_c = (double *)malloc(size);
		g_lambda = (double *)malloc(size);

		Dr_g_a = (double *)malloc(size);
		Dr_g_b = (double *)malloc(size);
		Dr_g_h = (double *)malloc(size);
		Dr_g_c = (double *)malloc(size);

		Dz_g_a = (double *)malloc(size);
		Dz_g_b = (double *)malloc(size);
		Dz_g_h = (double *)malloc(size);
		Dz_g_c = (double *)malloc(size);

		// DETERMINANT.
		hdet = (double *)malloc(size);
		Dr_hdet = (double *)malloc(size);
		Dz_hdet = (double *)malloc(size);

		// TRACELESS EXTRINSIC CURVATURE.
		A_a = (double *)malloc(size);
		A_b = (double *)malloc(size);
		A_h = (double *)malloc(size);
		A_c = (double *)malloc(size);
		A_lambda = (double *)malloc(size);

		Dr_A_a = (double *)malloc(size);
		Dr_A_b = (double *)malloc(size);
		Dr_A_h = (double *)malloc(size);
		Dr_A_c = (double *)malloc(size);
		Dr_A_lambda = (double *)malloc(size);

		Dz_A_a = (double *)malloc(size);
		Dz_A_b = (double *)malloc(size);
		Dz_A_h = (double *)malloc(size);
		Dz_A_c = (double *)malloc(size);
		Dz_A_lambda = (double *)malloc(size);

		// CONFORMAL FACTOR.
		phi = (double *)malloc(size);
		psi = (double *)malloc(size);
		psi4 = (double *)malloc(size);
		Dr_phi = (double *)malloc(size);
		Dz_phi = (double *)malloc(size);
		Drr_phi = (double *)malloc(size);
		Dzz_phi = (double *)malloc(size);
		Drz_phi = (double *)malloc(size);

		// TRACE OF EXTRINSIC CURVATURE.
		K = (double *)malloc(size);
		Dr_K = (double *)malloc(size);
		Dz_K = (double *)malloc(size);

		// DELTAS.
		Deltar = (double *)malloc(size);
		Dr_Deltar = (double *)malloc(size);
		Dz_Deltar = (double *)malloc(size);

		Deltaz = (double *)malloc(size);
		Dr_Deltaz = (double *)malloc(size);
		Dz_Deltaz = (double *)malloc(size);

		// REGULARIZATION.
		DDeltar_r = (double *)malloc(size);
		DDalpha_r = (double *)malloc(size);
		DDphi_r = (double *)malloc(size);

		// RICCI TENSOR.
		R_a = (double *)malloc(size);
		R_b = (double *)malloc(size);
		R_h = (double *)malloc(size);
		R_c = (double *)malloc(size);
		R_lambda = (double *)malloc(size);
		RSCAL = (double *)malloc(size);

		// SECOND COVARIANT DERIVATIVE OF THE LAPSE
		D2alpha_a = (double *)malloc(size);
		D2alpha_b = (double *)malloc(size);
		D2alpha_h = (double *)malloc(size);
		D2alpha_c = (double *)malloc(size);
		D2alpha_lambda = (double *)malloc(size);
		LaplaAlpha = (double *)malloc(size);

		// QUADRATIC QUANTITIES OF A.
		A2_a = (double *)malloc(size);
		A2_b = (double *)malloc(size);
		A2_h = (double *)malloc(size);
		A2_c = (double *)malloc(size);
		A2_lambda = (double *)malloc(size);
		A2 = (double *)malloc(size);

		// CONTRAVARIANT A.
		upA_a = (double *)malloc(size);
		upA_b = (double *)malloc(size);
		upA_h = (double *)malloc(size);
		upA_c = (double *)malloc(size);

		// DIVERGENCE OF A.
		divA_r = (double *)malloc(size);
		divA_z = (double *)malloc(size);

		// TERM IN DELTA SOURCES.
		ADelta_r = (double *)malloc(size);
		ADelta_z = (double *)malloc(size);

		// CONSTRAINTS.
		ham = (double *)malloc(size);
		mom_r = (double *)malloc(size);
		mom_z = (double *)malloc(size);
		CDeltar = (double *)malloc(size);
		CDeltaz = (double *)malloc(size);
		Clambda = (double *)malloc(size);
		CA_lambda = (double *)malloc(size);

		// OLD ARRAYS.
		alpha_p = (double *)malloc(size);
		phi_p = (double *)malloc(size);
		a_p = (double *)malloc(size);
		b_p = (double *)malloc(size);
		h_p = (double *)malloc(size);
		c_p = (double *)malloc(size);
		lambda_p = (double *)malloc(size);
		K_p = (double *)malloc(size);
		A_a_p = (double *)malloc(size);
		A_b_p = (double *)malloc(size);
		A_h_p = (double *)malloc(size);
		A_c_p = (double *)malloc(size);
		A_lambda_p = (double *)malloc(size);
		Deltar_p = (double *)malloc(size);
		Deltaz_p = (double *)malloc(size);

		// ACCUMULATION ARRAYS.
		alpha_a = (double *)malloc(size);
		phi_a = (double *)malloc(size);
		a_a = (double *)malloc(size);
		b_a = (double *)malloc(size);
		h_a = (double *)malloc(size);
		c_a = (double *)malloc(size);
		lambda_a = (double *)malloc(size);
		K_a = (double *)malloc(size);
		A_a_a = (double *)malloc(size);
		A_b_a = (double *)malloc(size);
		A_h_a = (double *)malloc(size);
		A_c_a = (double *)malloc(size);
		A_lambda_a = (double *)malloc(size);
		Deltar_a = (double *)malloc(size);
		Deltaz_a = (double *)malloc(size);

		// SOURCES.
		salpha = (double *)malloc(size);
		sphi = (double *)malloc(size);
		sa = (double *)malloc(size);
		sb = (double *)malloc(size);
		sh = (double *)malloc(size);
		sc = (double *)malloc(size);
		slambda = (double *)malloc(size);
		sK = (double *)malloc(size);
		sA_a = (double *)malloc(size);
		sA_b = (double *)malloc(size);
		sA_h = (double *)malloc(size);
		sA_c = (double *)malloc(size);
		sA_lambda = (double *)malloc(size);
		sDeltar = (double *)malloc(size);
		sDeltaz = (double *)malloc(size);

		// ELLIPTIC COEFFICENTS
		ell_a = (double *)malloc(size);
		ell_b = (double *)malloc(size);
		ell_c = (double *)malloc(size);
		ell_d = (double *)malloc(size);
		ell_e = (double *)malloc(size);
		ell_s = (double *)malloc(size);
		ell_f = (double *)malloc(size);

		// MIXED DERIVATIVES FOR BICUBIC INTERPOLATION.
		Drz_A_a = (double *)malloc(size);
		Drz_A_b = (double *)malloc(size);
		Drz_A_h = (double *)malloc(size);
		Drz_A_c = (double *)malloc(size);
		Drz_K = (double *)malloc(size);
		Drrz_a = (double *)malloc(size);
		Drrz_b = (double *)malloc(size);
		Drrz_h = (double *)malloc(size);
		Drrz_c = (double *)malloc(size);
		Drrz_phi = (double *)malloc(size);
		Drzz_a = (double *)malloc(size);
		Drzz_b = (double *)malloc(size);
		Drzz_h = (double *)malloc(size);
		Drzz_c = (double *)malloc(size);
		Drzz_phi = (double *)malloc(size);

		// PHYSICAL METRIC.
		phys_a = (double *)malloc(size);
		phys_b = (double *)malloc(size);
		phys_h = (double *)malloc(size);
		phys_c = (double *)malloc(size);
		phys_lambda = (double *)malloc(size);

		// PHYSICAL SPHERICAL METRIC.
		g_rr = (double *)malloc(size);
		g_rth = (double *)malloc(size);
		g_thth = (double *)malloc(size);
		g_phph = (double *)malloc(size);
		ig_rr = (double *)malloc(size);
		ig_rth = (double *)malloc(size);
		ig_thth = (double *)malloc(size);
		ig_phph = (double *)malloc(size);

		// SPEED OF LIGHT.
		c_rr = (double *)malloc(size);
		c_th = (double *)malloc(size);
		c_ph = (double *)malloc(size);

		// COURANT CONDITION.
		courant_rr = (double *)malloc(size);
		courant_th = (double *)malloc(size);
		courant_ph = (double *)malloc(size);

		// SHIFT.
		beta_r = (double *)malloc(size);
		beta_z = (double *)malloc(size);
		dtbeta_r = (double *)malloc(size);
		dtbeta_z = (double *)malloc(size);
		Dr_beta_r = (double *)malloc(size);
		Dz_beta_r = (double *)malloc(size);
		Drr_beta_r = (double *)malloc(size);
		Dzz_beta_r = (double *)malloc(size);
		Drz_beta_r = (double *)malloc(size);
		Dr_beta_z = (double *)malloc(size);
		Dz_beta_z = (double *)malloc(size);
		Drr_beta_z = (double *)malloc(size);
		Dzz_beta_z = (double *)malloc(size);
		Drz_beta_z = (double *)malloc(size);
		DD_beta_rr = (double *)malloc(size);
		Dr_dtbeta_r = (double *)malloc(size);
		Dz_dtbeta_r = (double *)malloc(size);
		Dr_dtbeta_z = (double *)malloc(size);
		Dz_dtbeta_z = (double *)malloc(size);
		div_beta = (double *)malloc(size);
		Dr_div_beta = (double *)malloc(size);
		Dz_div_beta = (double *)malloc(size);
		sbeta_r = (double *)malloc(size);
		sbeta_z = (double *)malloc(size);
		sdtbeta_r = (double *)malloc(size);
		sdtbeta_z = (double *)malloc(size);
		beta_r_p = (double *)malloc(size);
		beta_z_p = (double *)malloc(size);
		dtbeta_r_p = (double *)malloc(size);
		dtbeta_z_p = (double *)malloc(size);
		beta_r_a = (double *)malloc(size);
		beta_z_a = (double *)malloc(size);
		dtbeta_r_a = (double *)malloc(size);
		dtbeta_z_a = (double *)malloc(size);

		// ADVECTIVE DERIVATIVES.
		DAr_alpha = (double *)malloc(size);
		DAz_alpha = (double *)malloc(size);
		DAr_phi = (double *)malloc(size);
		DAz_phi = (double *)malloc(size);
		DAr_a = (double *)malloc(size);
		DAz_a = (double *)malloc(size);
		DAr_b = (double *)malloc(size);
		DAz_b = (double *)malloc(size);
		DAr_h = (double *)malloc(size);
		DAz_h = (double *)malloc(size);
		DAr_c = (double *)malloc(size);
		DAz_c = (double *)malloc(size);
		DAr_lambda = (double *)malloc(size);
		DAz_lambda = (double *)malloc(size);
		DAr_A_lambda = (double *)malloc(size);
		DAz_A_lambda = (double *)malloc(size);
		DAr_beta_r = (double *)malloc(size);
		DAz_beta_r = (double *)malloc(size);
		DAr_beta_z = (double *)malloc(size);
		DAz_beta_z = (double *)malloc(size);
		DAr_dtbeta_r = (double *)malloc(size);
		DAz_dtbeta_r = (double *)malloc(size);
		DAr_dtbeta_z = (double *)malloc(size);
		DAz_dtbeta_z = (double *)malloc(size);
		DAr_K = (double *)malloc(size);
		DAz_K = (double *)malloc(size);
		DAr_A_a = (double *)malloc(size);
		DAz_A_a = (double *)malloc(size);
		DAr_A_b = (double *)malloc(size);
		DAz_A_b = (double *)malloc(size);
		DAr_A_h = (double *)malloc(size);
		DAz_A_h = (double *)malloc(size);
		DAr_A_c = (double *)malloc(size);
		DAz_A_c = (double *)malloc(size);
		DAr_Deltar = (double *)malloc(size);
		DAz_Deltar = (double *)malloc(size);
		DAr_Deltaz = (double *)malloc(size);
		DAz_Deltaz = (double *)malloc(size);
	}
	// Deallocate arrays.
	else
	{
		// AUXILIARY ARRAYS.
		free(auxarray);
		free(auxarray2);
		free(auxarray3);
		free(auxarray4);

		// COORDINATE GRIDS.
		free(r);
		free(z);
		free(rr);
		free(r2);

		// LASPE.
		free(alpha);
		free(Dr_alpha);
		free(Dz_alpha);
		free(Drr_alpha);
		free(Dzz_alpha);
		free(Drz_alpha);

		// BONA-MASSO SOURCE.
		free(falpha);

		// CONFORMAL METRIC.
		free(a);
		free(b);
		free(h);
		free(c);
		free(lambda);

		free(Dr_a);
		free(Dr_b);
		free(Dr_h);
		free(Dr_c);
		free(Dr_lambda);

		free(Dz_a);
		free(Dz_b);
		free(Dz_h);
		free(Dz_c);
		free(Dz_lambda);

		free(Drr_a);
		free(Drr_b);
		free(Drr_h);
		free(Drr_c);
		free(Drr_lambda);

		free(Dzz_a);
		free(Dzz_b);
		free(Dzz_h);
		free(Dzz_c);
		free(Dzz_lambda);

		free(Drz_a);
		free(Drz_b);
		free(Drz_h);
		free(Drz_c);
		free(Drz_lambda);

		// INVERSE CONORMAL METRIC.
		free(g_a);
		free(g_b);
		free(g_h);
		free(g_c);
		free(g_lambda);

		free(Dr_g_a);
		free(Dr_g_b);
		free(Dr_g_h);
		free(Dr_g_c);

		free(Dz_g_a);
		free(Dz_g_b);
		free(Dz_g_h);
		free(Dz_g_c);

		// DETERMINANT.
		free(hdet);
		free(Dr_hdet);
		free(Dz_hdet);

		// TRACELESS EXTRINSIC CURVATURE.
		free(A_a);
		free(A_b);
		free(A_h);
		free(A_c);
		free(A_lambda);

		free(Dr_A_a);
		free(Dr_A_b);
		free(Dr_A_h);
		free(Dr_A_c);
		free(Dr_A_lambda);

		free(Dz_A_a);
		free(Dz_A_b);
		free(Dz_A_h);
		free(Dz_A_c);
		free(Dz_A_lambda);

		// CONFORMAL FACTOR.
		free(phi);
		free(psi);
		free(psi4);
		free(Dr_phi);
		free(Dz_phi);
		free(Drr_phi);
		free(Dzz_phi);
		free(Drz_phi);

		// TRACE OF EXTRINSIC CURVATURE.
		free(K);
		free(Dr_K);
		free(Dz_K);

		// DELTAS.
		free(Deltar);
		free(Dr_Deltar);
		free(Dz_Deltar);

		free(Deltaz);
		free(Dr_Deltaz);
		free(Dz_Deltaz);

		// REGULARIZATION.
		free(DDeltar_r);
		free(DDalpha_r);
		free(DDphi_r);

		// RICCI TENSOR.
		free(R_a);
		free(R_b);
		free(R_h);
		free(R_c);
		free(R_lambda);
		free(RSCAL);

		// SECOND COVARIANT DERIVATIVE OF THE LAPSE
		free(D2alpha_a);
		free(D2alpha_b);
		free(D2alpha_h);
		free(D2alpha_c);
		free(D2alpha_lambda);
		free(LaplaAlpha);

		// QUADRATIC QUANTITIES OF A.
		free(A2_a);
		free(A2_b);
		free(A2_h);
		free(A2_c);
		free(A2_lambda);
		free(A2);

		// CONTRAVARIANT A.
		free(upA_a);
		free(upA_b);
		free(upA_h);
		free(upA_c);

		// DIVERGENCE OF A.
		free(divA_r);
		free(divA_z);

		// TERM IN DELTA SOURCES.
		free(ADelta_r);
		free(ADelta_z);

		// CONSTRAINTS.
		free(ham);
		free(mom_r);
		free(mom_z);
		free(CDeltar);
		free(CDeltaz);
		free(Clambda);
		free(CA_lambda);

		// OLD ARRAYS.
		free(alpha_p);
		free(phi_p);
		free(a_p);
		free(b_p);
		free(h_p);
		free(c_p);
		free(lambda_p);
		free(K_p);
		free(A_a_p);
		free(A_b_p);
		free(A_h_p);
		free(A_c_p);
		free(A_lambda_p);
		free(Deltar_p);
		free(Deltaz_p);

		// ACCUMULATION ARRAYS.
		free(alpha_a);
		free(phi_a);
		free(a_a);
		free(b_a);
		free(h_a);
		free(c_a);
		free(lambda_a);
		free(K_a);
		free(A_a_a);
		free(A_b_a);
		free(A_h_a);
		free(A_c_a);
		free(A_lambda_a);
		free(Deltar_a);
		free(Deltaz_a);

		// SOURCES.
		free(salpha);
		free(sphi);
		free(sa);
		free(sb);
		free(sh);
		free(sc);
		free(slambda);
		free(sK);
		free(sA_a);
		free(sA_b);
		free(sA_h);
		free(sA_c);
		free(sA_lambda);
		free(sDeltar);
		free(sDeltaz);

		// ELLIPTIC COEFFICENTS
		free(ell_a);
		free(ell_b);
		free(ell_c);
		free(ell_d);
		free(ell_e);
		free(ell_s);
		free(ell_f);

		// MIXED DERIVATIVES FOR BICUBIC INTERPOLATION.
		free(Drz_A_a);
		free(Drz_A_b);
		free(Drz_A_h);
		free(Drz_A_c);
		free(Drz_K);
		free(Drrz_a);
		free(Drrz_b);
		free(Drrz_h);
		free(Drrz_c);
		free(Drrz_phi);
		free(Drzz_a);
		free(Drzz_b);
		free(Drzz_h);
		free(Drzz_c);
		free(Drzz_phi);

		// PHYSICAL METRIC.
		free(phys_a);
		free(phys_b);
		free(phys_h);
		free(phys_c);
		free(phys_lambda);

		// PHYSICAL SPHERICAL METRIC.
		free(g_rr);
		free(g_rth);
		free(g_thth);
		free(g_phph);
		free(ig_rr);
		free(ig_rth);
		free(ig_thth);
		free(ig_phph);

		// SPEED OF LIGHT.
		free(c_rr);
		free(c_th);
		free(c_ph);

		// COURANT CONDITION.
		free(courant_rr);
		free(courant_th);
		free(courant_ph);

		// SHIFT.
		free(beta_r);
		free(beta_z);
		free(dtbeta_r);
		free(dtbeta_z);
		free(Dr_beta_r);
		free(Dz_beta_r);
		free(Drr_beta_r);
		free(Dzz_beta_r);
		free(Drz_beta_r);
		free(Dr_beta_z);
		free(Dz_beta_z);
		free(Drr_beta_z);
		free(Dzz_beta_z);
		free(Drz_beta_z);
		free(DD_beta_rr);
		free(Dr_dtbeta_r);
		free(Dz_dtbeta_r);
		free(Dr_dtbeta_z);
		free(Dz_dtbeta_z);
		free(div_beta);
		free(Dr_div_beta);
		free(Dz_div_beta);
		free(sbeta_r);
		free(sbeta_z);
		free(sdtbeta_r);
		free(sdtbeta_z);
		free(beta_r_p);
		free(beta_z_p);
		free(dtbeta_r_p);
		free(dtbeta_z_p);
		free(beta_r_a);
		free(beta_z_a);
		free(dtbeta_r_a);
		free(dtbeta_z_a);

		// ADVECTIVE DERIVATIVES.
		free(DAr_alpha);
		free(DAz_alpha);
		free(DAr_phi);
		free(DAz_phi);
		free(DAr_a);
		free(DAz_a);
		free(DAr_b);
		free(DAz_b);
		free(DAr_h);
		free(DAz_h);
		free(DAr_c);
		free(DAz_c);
		free(DAr_lambda);
		free(DAz_lambda);
		free(DAr_A_lambda);
		free(DAz_A_lambda);
		free(DAr_beta_r);
		free(DAz_beta_r);
		free(DAr_beta_z);
		free(DAz_beta_z);
		free(DAr_dtbeta_r);
		free(DAz_dtbeta_r);
		free(DAr_dtbeta_z);
		free(DAz_dtbeta_z);
		free(DAr_K);
		free(DAz_K);
		free(DAr_A_a);
		free(DAz_A_a);
		free(DAr_A_b);
		free(DAz_A_b);
		free(DAr_A_h);
		free(DAz_A_h);
		free(DAr_A_c);
		free(DAz_A_c);
		free(DAr_Deltar);
		free(DAz_Deltar);
		free(DAr_Deltaz);
		free(DAz_Deltaz);
	}
}
