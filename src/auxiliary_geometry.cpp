#include "tools.h"

#include "param.h"
#include "arrays.h"

#include "geometry_A2.h"
#include "geometry_ADelta.h"
#include "geometry_alpha_derivatives.h"
#include "geometry_alpha_second_covariant_derivative.h"
#include "geometry_detrace.h"
#include "geometry_divA.h"
#include "geometry_inverse_metric.h"
#include "geometry_ricci.h"
#include "geometry_upA.h"

#include "derivatives.h"

#include "shift_divergence.h"

#undef DEBUG

void auxiliary_geometry(void)
{
	// Auxiliary integers.
	int i, j, k;

	// Preilimaries:
	// First calculate inverse metric and g_lambda.
	#pragma omp parallel shared(g_a, g_b, g_h, g_c, g_lambda, hdet)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			// Inverse metric.
			calculate_inverse_metric(g_a + k, g_b + k, g_h + k, g_c + k,
				g_lambda + k, hdet + k,
				r[k], a[k], b[k], h[k], c[k], lambda[k]);
		}
	}

	// Calculate conformal factor, but only if not initial time,
	// and not maximal slicing, where it does not evolve.
	//if ((t > 0.) && (!(strcmp(slicing, "maximal") == 0) || !(strcmp(shift, "none") == 0)))
	#pragma omp parallel shared(psi, psi4)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			psi[k] = exp(phi[k]);
			psi4[k] = psi[k] * psi[k] * psi[k] * psi[k];
		}
	}

	// Conformal factor derivatives.
	diff1r(Dr_phi, phi, 1);
	diff1z(Dz_phi, phi, 1);

	// Regularization.
	#pragma omp parallel shared(auxarray)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			auxarray[k] = Dr_phi[k] / r[k];
		}
	}
	diff1r(DDphi_r, auxarray, 1);

	// Second derivatives.
	diff2r(Drr_phi, phi, 1);
	diff2z(Dzz_phi, phi, 1);
	diff2rz(Drz_phi, phi, 1, 1);

	// Calculate derivatives of K if not maximal slicing.
	if (!(strcmp(slicing, "maximal") == 0))
	{
		diff1r(Dr_K, K, 1);
		diff1z(Dz_K, K, 1);
	}

	// Derivatives of the conformal metric.
	diff1r(Dr_a, a, 1);
	diff1r(Dr_b, b, 1);
	diff1r(Dr_h, h, 1);
	diff1r(Dr_c, c, 1);
	diff1r(Dr_lambda, lambda, 1);
	diff1z(Dz_a, a, 1);
	diff1z(Dz_b, b, 1);
	diff1z(Dz_h, h, 1);
	diff1z(Dz_c, c, -1);
	diff1z(Dz_lambda, lambda, 1);

	// Derivatives of the inverse conformal metric.
	diff1r(Dr_g_a, g_a, 1);
	diff1r(Dr_g_b, g_b, 1);
	diff1r(Dr_g_h, g_h, 1);
	diff1r(Dr_g_c, g_c, 1);
	diff1z(Dz_g_a, g_a, 1);
	diff1z(Dz_g_b, g_b, 1);
	diff1z(Dz_g_h, g_h, 1);
	diff1z(Dz_g_c, g_c, -1);

	// Derivatives of the traceless extrinsic curvature.
	diff1r(Dr_A_a, A_a, 1);
	diff1r(Dr_A_b, A_b, 1);
	diff1r(Dr_A_h, A_h, 1);
	diff1r(Dr_A_c, A_c, 1);
	diff1r(Dr_A_lambda, A_lambda, 1);
	diff1z(Dz_A_a, A_a, 1);
	diff1z(Dz_A_b, A_b, 1);
	diff1z(Dz_A_h, A_h, 1);
	diff1z(Dz_A_c, A_c, -1);
	diff1z(Dz_A_lambda, A_lambda, 1);

	// Delta derivatives.
	diff1r(Dr_Deltar, Deltar, -1);
	diff1r(Dr_Deltaz, Deltaz, 1);
	diff1z(Dz_Deltar, Deltar, 1);
	diff1z(Dz_Deltaz, Deltaz, -1);
	// Regularization.
	#pragma omp parallel shared(auxarray)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			auxarray[k] = Deltar[k] / r[k];
		}
	}
	diff1r(DDeltar_r, auxarray, 1);

	// Calculate second derivatives.
	// Conformal metric.
	diff2r(Drr_a, a, 1);
	diff2z(Dzz_a, a, 1);
	diff2rz(Drz_a, a, 1, 1);
	diff2r(Drr_b, b, 1);
	diff2z(Dzz_b, b, 1);
	diff2rz(Drz_b, b, 1, 1);
	diff2r(Drr_h, h, 1);
	diff2z(Dzz_h, h, 1);
	diff2rz(Drz_h, h, 1, 1);
	diff2r(Drr_c, c, 1);
	diff2z(Dzz_c, c, -1);
	diff2rz(Drz_c, c, 1, -1);
	diff2r(Drr_lambda, lambda, 1);
	diff2z(Dzz_lambda, lambda, 1);
	diff2rz(Drz_lambda, lambda, 1, 1);

	// Conformal Ricci.
	#pragma omp parallel shared(R_a, R_b, R_h, R_c, R_lambda, RSCAL)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			// Calculate Conformal Ricci tensor.
			calculate_ricci(R_a + k, R_b + k, R_h + k, R_c + k, R_lambda + k, RSCAL + k,
				r[k], a[k], b[k], h[k], c[k], lambda[k],
				g_a[k], g_b[k], g_h[k], g_c[k], g_lambda[k],
				Dr_a[k], Dr_b[k], Dr_h[k], Dr_c[k], Dr_lambda[k],
				Dz_a[k], Dz_b[k], Dz_h[k], Dz_c[k], Dz_lambda[k],
				Drr_a[k], Drr_b[k], Drr_h[k], Drr_c[k], Drr_lambda[k],
				Dzz_a[k], Dzz_b[k], Dzz_h[k], Dzz_c[k], Dzz_lambda[k],
				Drz_a[k], Drz_b[k], Drz_h[k], Drz_c[k], Drz_lambda[k],
				Deltar[k], Deltaz[k],
				Dr_Deltar[k], Dr_Deltaz[k], Dz_Deltar[k], Dz_Deltaz[k],
				DDeltar_r[k]);
		}
	}

#ifdef DEBUG
	writeSingleFile(R_a, "conf_R_a.asc");
	writeSingleFile(R_b, "conf_R_b.asc");
	writeSingleFile(R_h, "conf_R_h.asc");
	writeSingleFile(R_c, "conf_R_c.asc");
	writeSingleFile(R_lambda, "conf_R_lambda.asc");
	writeSingleFile(RSCAL, "conf_RSCAL.asc");
#endif

	// Detrace the traceless extrinsic curvature of any numerical 
	// error that may have accumulated.
	detrace_A();

	// Main calculations.
	#pragma omp parallel shared(R_a, R_b, R_h, R_c, R_lambda, RSCAL,\
	A2_a, A2_b, A2_h, A2_c, A2_lambda, A2,\
	upA_a, upA_b, upA_h, upA_c, divA_r, divA_z, ADelta_r, ADelta_z)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			// Calculate physical Ricci tensor.
			calculate_phys_ricci(R_a + k, R_b + k, R_h + k, R_c + k,
				R_lambda + k, RSCAL + k,
				r[k], a[k], b[k], h[k], c[k], lambda[k],
				g_a[k], g_b[k], g_h[k], g_c[k], g_lambda[k],
				Dr_a[k], Dr_b[k], Dr_h[k], Dr_c[k], Dr_lambda[k],
				Dz_a[k], Dz_b[k], Dz_h[k], Dz_c[k], Dz_lambda[k],
				Dr_phi[k], Dz_phi[k], Drr_phi[k], Dzz_phi[k], Drz_phi[k],
				DDphi_r[k], psi4[k]);

			// Calculate quadratic quantities.
			calculate_A2(A2_a + k, A2_b + k, A2_h + k, A2_c + k, A2_lambda + k, A2 + k,
				r[k], A_a[k], A_b[k], A_h[k], A_c[k], A_lambda[k],
				g_a[k], g_b[k], g_h[k], g_c[k], g_lambda[k]);

			// Calculate contravariant traceless extrinsic curvature.
			calculate_upA(upA_a + k, upA_b + k, upA_h + k, upA_c + k,
				r[k], g_a[k], g_b[k], g_h[k], g_c[k],
				A_a[k], A_b[k], A_h[k], A_c[k]);

			// Calculate divergence.
			calculate_divA(divA_r + k, divA_z + k,
				r[k], a[k], b[k], h[k], c[k],
				A_a[k], A_b[k], A_h[k], A_c[k], A_lambda[k],
				g_a[k], g_b[k], g_h[k], g_c[k], g_lambda[k],
				Dr_a[k], Dr_b[k], Dr_h[k], Dr_c[k],
				Dz_a[k], Dz_b[k], Dz_h[k], Dz_c[k],
				Dr_A_a[k], Dr_A_b[k], Dr_A_h[k], Dr_A_c[k],
				Dz_A_a[k], Dz_A_b[k], Dz_A_h[k], Dz_A_c[k]);

			// Calculate term in Delta sources.
			calculate_ADelta(ADelta_r + k, ADelta_z + k,
				r[k], a[k], b[k], h[k], c[k], lambda[k],
				g_a[k], g_b[k], g_h[k], g_c[k],
				A_a[k], A_b[k], A_h[k], A_c[k],
				Dr_a[k], Dr_b[k], Dr_h[k], Dr_c[k],
				Dz_a[k], Dz_b[k], Dz_h[k], Dz_c[k]);
		}
	}

	return;
}