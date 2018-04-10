void calculate_sources_geometry_maximal(double *sa, double *sb, double *sh, double *sc, double *slambda,
	double *sA_a, double *sA_b, double *sA_h, double *sA_c, double *sA_lambda, double *sphi, double *sK,
	double *sDeltar, double *sDeltaz,
	const double r, const double alpha, const double A_a, const double A_b, const double A_h, const double A_c, const double A_lambda,
	const double D2alpha_a, const double D2alpha_b, const double D2alpha_h, const double D2alpha_c, const double D2alpha_lambda,
	const double R_a, const double R_b, const double R_h, const double R_c, const double R_lambda, const double psi4,
	const double A2_a, const double A2_b, const double A2_h, const double A2_c, const double A2_lambda,
	const double upA_a, const double upA_b, const double upA_c, const double Dr_alpha, const double Dz_alpha,
	const double divA_r, const double divA_z, const double ADelta_r, const double ADelta_z,
	const double Dr_phi, const double Dz_phi, const double eta)
{
	const double third = 1. / 3.;

	// Sources for the conformal metric.
	*sa = -2. * alpha * A_a;
	*sb = -2. * alpha * A_b;
	*sh = -2. * alpha * A_h;
	*sc = -2. * alpha * A_c;
	*slambda = -2. * alpha * A_lambda;

	// We can simplify the expressions if maximal slicing has been chosen.
	// Since K = Dr_K = Dz_K = Dt_K = 0 and LaplaAlpha = alpha R.
	// Sources for the traceless extrinsic curvature.
	*sA_a = (-D2alpha_a + alpha * R_a) / psi4
		+ alpha * (-2. * A2_a);
	*sA_b = (-D2alpha_b + alpha * R_b) / psi4
		+ alpha * (-2. * A2_b);
	*sA_h = (-D2alpha_h + alpha * R_h) / psi4
		+ alpha * (-2. * A2_h);
	*sA_c = (-D2alpha_c + alpha * R_c) / psi4
		+ alpha * (-2. * A2_c);
	*sA_lambda = (-D2alpha_lambda + alpha * R_lambda) / psi4
		+ alpha * (-2. * A2_lambda);

	// Sources for conformal factor.
	*sphi = 0.;

	// Sources for trace of the extrinsic curvature.
	*sK = 0.;

	// Sources for Deltas.
	*sDeltar = -2. * (upA_a * Dr_alpha + r * upA_c * Dz_alpha)
		- alpha * (2. - eta) * divA_r
		+ 2. * alpha * ADelta_r
		+ alpha * eta * 6. * (upA_a * Dr_phi + r * upA_c * Dz_phi);
	*sDeltaz = -2. * (upA_b * Dz_alpha + r * upA_c * Dr_alpha)
		- alpha * (2. - eta) * divA_z
		+ 2. * alpha * ADelta_z
		+ alpha * eta * 6. * (upA_b * Dz_phi + r * upA_c * Dr_phi);
}

void calculate_sources_geometry_full(double *sa, double *sb, double *sh, double *sc, double *slambda,
	double *sA_a, double *sA_b, double *sA_h, double *sA_c, double *sA_lambda, double *sphi, double *sK,
	double *sDeltar, double *sDeltaz,
	const double r, const double alpha, const double a, const double b, const double h, const double c, const double lambda,
	const double A_a, const double A_b, const double A_h, const double A_c, const double A_lambda,
	const double D2alpha_a, const double D2alpha_b, const double D2alpha_h, const double D2alpha_c, const double D2alpha_lambda, const double LaplaAlpha,
	const double R_a, const double R_b, const double R_h, const double R_c, const double R_lambda, const double R,
	const double A2_a, const double A2_b, const double A2_h, const double A2_c, const double A2_lambda, const double A2,
	const double upA_a, const double upA_b, const double upA_c,
	const double ADelta_r, const double ADelta_z, const double Dr_alpha, const double Dz_alpha,
	const double Dr_phi, const double Dz_phi, const double psi4, const double divA_r, const double divA_z,
	const double K, const double Dr_K, const double Dz_K,
	const double g_a, const double g_b, const double g_c, const double eta)
{
	const double third = 1. / 3.;
	const double sixth = 1. / 6.;

	// Sources for the conformal metric.
	*sa = -2. * alpha * A_a;
	*sb = -2. * alpha * A_b;
	*sh = -2. * alpha * A_h;
	*sc = -2. * alpha * A_c;
	*slambda = -2. * alpha * A_lambda;

	// Sources for the traceless extrinsic curvature.
	*sA_a = (-D2alpha_a + alpha * R_a) / psi4 - third * a * (-LaplaAlpha + alpha * R)
		+ alpha * (K * A_a - 2. * A2_a);
	*sA_b = (-D2alpha_b + alpha * R_b) / psi4 - third * b * (-LaplaAlpha + alpha * R)
		+ alpha * (K * A_b - 2. * A2_b);
	*sA_h = (-D2alpha_h + alpha * R_h) / psi4 - third * h * (-LaplaAlpha + alpha * R)
		+ alpha * (K * A_h - 2. * A2_h);
	*sA_c = (-D2alpha_c + alpha * R_c) / psi4 - third * c * (-LaplaAlpha + alpha * R)
		+ alpha * (K * A_c - 2. * A2_c);
	*sA_lambda = (-D2alpha_lambda + alpha * R_lambda) / psi4 - third * lambda * (-LaplaAlpha + alpha * R)
		+ alpha * (K * A_lambda - 2. * A2_lambda);

	// Sources for conformal factor.
	*sphi = -sixth * alpha * K;

	// Sources for trace of the extrinsic curvature.
	*sK = -LaplaAlpha + alpha * (A2 + third * K * K);

	// Sources for Deltas.
	*sDeltar = -2. * (upA_a * Dr_alpha + r * upA_c * Dz_alpha)
		- alpha * (2. - eta) * divA_r
		+ 2. * alpha * ADelta_r
		+ alpha * eta * (6. * (upA_a * Dr_phi + r * upA_c * Dz_phi) - 2. * third * (g_a * Dr_K + r * g_c * Dz_K));
	*sDeltaz = -2. * (upA_b * Dz_alpha + r * upA_c * Dr_alpha)
		- alpha * (2. - eta) * divA_z
		+ 2. * alpha * ADelta_z
		+ alpha * eta * (6. * (upA_b * Dz_phi + r * upA_c * Dr_phi) - 2. * third * (g_b * Dz_K + r * g_c * Dr_K));
}

#include "tools.h"

#include "param.h"
#include "arrays.h"

#define FULL

void calculate_sources_geometry_shift(double *sa, double *sb, double *sh, double *sc, double *slambda,
	double *sA_a, double *sA_b, double *sA_h, double *sA_c, double *sA_lambda, double *sphi, double *sK,
	double *sDeltar, double *sDeltaz,
	const double r, const double a, const double b, const double h, const double c, const double lambda,
	const double A_a, const double A_b, const double A_h, const double A_c, const double A_lambda,
	const double beta_r, const double Dr_beta_r, const double Dz_beta_r, const double Drr_beta_r, const double Dzz_beta_r, const double Drz_beta_r,
	const double beta_z, const double Dr_beta_z, const double Dz_beta_z, const double Drr_beta_z, const double Dzz_beta_z, const double Drz_beta_z,
	const double DD_beta_rr, const double div_beta, const double Dr_div_beta, const double Dz_div_beta,
	const double DAr_a, const double DAr_b, const double DAr_h, const double DAr_c, const double DAr_lambda,
	const double DAz_a, const double DAz_b, const double DAz_h, const double DAz_c, const double DAz_lambda,
	const double DAr_A_a, const double DAr_A_b, const double DAr_A_h, const double DAr_A_c, const double DAr_A_lambda,
	const double DAz_A_a, const double DAz_A_b, const double DAz_A_h, const double DAz_A_c, const double DAz_A_lambda,
	const double DAr_phi, const double DAz_phi, const double DAr_K, const double DAz_K,
	const double g_a, const double g_b, const double g_h, const double g_c,
	const double Deltar, const double DAr_Deltar, const double DAz_Deltar,
	const double Deltaz, const double DAr_Deltaz, const double DAz_Deltaz,
	const double sigma);

void sources_geometry(void)
{
	int k;

	if (strcmp(slicing, "maximal") == 0)
	{
		#pragma omp parallel shared(sa, sb, sh, sc, slambda,\
		sA_a, sA_b, sA_h, sA_c, sA_lambda, sphi, sK, sDeltar, sDeltaz)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
#ifdef FULL
				calculate_sources_geometry_full(sa + k, sb + k, sh + k, sc + k, slambda + k,
					sA_a + k, sA_b + k, sA_h + k, sA_c + k, sA_lambda + k,
					sphi + k, sK + k, sDeltar + k, sDeltaz + k,
					r[k], alpha[k], a[k], b[k], h[k], c[k], lambda[k],
					A_a[k], A_b[k], A_h[k], A_c[k], A_lambda[k],
					D2alpha_a[k], D2alpha_b[k], D2alpha_h[k], D2alpha_c[k], D2alpha_lambda[k], LaplaAlpha[k],
					R_a[k], R_b[k], R_h[k], R_c[k], R_lambda[k], RSCAL[k],
					A2_a[k], A2_b[k], A2_h[k], A2_c[k], A2_lambda[k], A2[k],
					upA_a[k], upA_b[k], upA_c[k], ADelta_r[k], ADelta_z[k],
					Dr_alpha[k], Dz_alpha[k], Dr_phi[k], Dz_phi[k], psi4[k],
					divA_r[k], divA_z[k], 0.0, 0.0, 0.0, g_a[k], g_b[k], g_c[k], eta);

				// Set to zero for conformal factor and curvature trace.
				sphi[k] = 0.0;
				sK[k] = 0.0;
#else
				calculate_sources_geometry_maximal(sa + k, sb + k, sh + k, sc + k, slambda + k,
					sA_a + k, sA_b + k, sA_h + k, sA_c + k, sA_lambda + k,
					sphi + k, sK + k, sDeltar + k, sDeltaz + k,
					r[k], alpha[k], A_a[k], A_b[k], A_h[k], A_c[k], A_lambda[k],
					D2alpha_a[k], D2alpha_b[k], D2alpha_h[k], D2alpha_c[k], D2alpha_lambda[k],
					R_a[k], R_b[k], R_h[k], R_c[k], R_lambda[k], psi4[k],
					A2_a[k], A2_b[k], A2_h[k], A2_c[k], A2_lambda[k],
					upA_a[k], upA_b[k], upA_c[k], Dr_alpha[k], Dz_alpha[k],
					divA_r[k], divA_z[k], ADelta_r[k], ADelta_z[k],
					Dr_phi[k], Dz_phi[k], eta);
#endif
			}
		}
	}
	else
	{
		#pragma omp parallel shared(sa, sb, sh, sc, slambda,\
		sA_a, sA_b, sA_h, sA_c, sA_lambda, sphi, sK, sDeltar, sDeltaz)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				calculate_sources_geometry_full(sa + k, sb + k, sh + k, sc + k, slambda + k,
					sA_a + k, sA_b + k, sA_h + k, sA_c + k, sA_lambda + k,
					sphi + k, sK + k, sDeltar + k, sDeltaz + k,
					r[k], alpha[k], a[k], b[k], h[k], c[k], lambda[k],
					A_a[k], A_b[k], A_h[k], A_c[k], A_lambda[k],
					D2alpha_a[k], D2alpha_b[k], D2alpha_h[k], D2alpha_c[k], D2alpha_lambda[k], LaplaAlpha[k],
					R_a[k], R_b[k], R_h[k], R_c[k], R_lambda[k], RSCAL[k],
					A2_a[k], A2_b[k], A2_h[k], A2_c[k], A2_lambda[k], A2[k],
					upA_a[k], upA_b[k], upA_c[k], ADelta_r[k], ADelta_z[k],
					Dr_alpha[k], Dz_alpha[k], Dr_phi[k], Dz_phi[k], psi4[k],
					divA_r[k], divA_z[k], K[k], Dr_K[k], Dz_K[k], 
					g_a[k], g_b[k], g_c[k], eta);
			}
		}
	}
}