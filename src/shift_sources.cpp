#define FT1 0.5
#define FT2 0.5

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
	const double sigma)
{
	const double third = 1.0 / 3.0;
	const double sixth = 1.0 / 6.0;

	//     m    ^     ^        m   ^        m              ^  ^     m
	// beta  d  g   + g  d beta  + g  d beta  - (2/3)sigma g  D beta
	//        m  ij    mi j         mj i                    ij m
	*sa += beta_r * DAr_a + beta_z * DAz_a + 2.0 * (a * Dr_beta_r + r * c * Dr_beta_z) - 2.0 * third * sigma * a * div_beta;
	*sb += beta_r * DAr_b + beta_z * DAz_b + 2.0 * (b * Dz_beta_z + r * c * Dz_beta_r) - 2.0 * third * sigma * b * div_beta;
	*sh += beta_r * DAr_h + beta_z * DAz_h + 2.0 * h * (beta_r / r) - 2.0 * third * sigma * h * div_beta;
	*sc += beta_r * DAr_c + c * (beta_r / r) + beta_z * DAz_c + c * (Dr_beta_r + Dz_beta_z) + b * (Dr_beta_z / r) + a * (Dz_beta_r / r) - 2.0 * third * sigma * c * div_beta;
	*slambda += beta_r * DAr_lambda + beta_z * DAz_lambda + 2.0 * lambda * (beta_r / r) + 2.0 * c * (Dr_beta_z / r) - 2.0 * third * sigma * lambda * div_beta
		// Regularization.
		+ FT1 * (2.0 * h * (DD_beta_rr / r) + 2.0 * lambda * Dr_beta_r)
		+ (1.0 - FT1) * (2.0 * a * (DD_beta_rr / r) + 2.0 * lambda * (beta_r / r));

	//     m                       ^     m
	// beta  d  phi + (1/6)sigma * D beta
	//        m                     m
	*sphi += beta_r * DAr_phi + beta_z * DAz_phi + sixth * sigma * div_beta;

	//     m
	// beta  d K
	//        m
	*sK += beta_r * DAr_K + beta_z * DAz_K;

	//     m    ^     ^        m   ^        m              ^  ^     m
	// beta  d  A   + A  d beta  + A  d beta  - (2/3)sigma A  D beta
	//        m  ij    mi j         mj i                    ij m
	*sA_a += beta_r * DAr_A_a + beta_z * DAz_A_a + 2.0 * (A_a * Dr_beta_r + r * A_c * Dr_beta_z) - 2.0 * third * sigma * A_a * div_beta;
	*sA_b += beta_r * DAr_A_b + beta_z * DAz_A_b + 2.0 * (A_b * Dz_beta_z + r * A_c * Dz_beta_r) - 2.0 * third * sigma * A_b * div_beta;
	*sA_h += beta_r * DAr_A_h + beta_z * DAz_A_h + 2.0 * A_h * (beta_r / r) - 2.0 * third * sigma * A_h * div_beta;
	*sA_c += beta_r * DAr_A_c + A_c * (beta_r / r) + beta_z * DAz_A_c + A_c * (Dr_beta_r + Dz_beta_z) + A_b * (Dr_beta_z / r) + A_a * (Dz_beta_r / r) - 2.0 * third * sigma * A_c * div_beta;
	*sA_lambda += beta_r * DAr_A_lambda + beta_z * DAz_A_lambda + 2.0 * A_lambda * (beta_r / r) + 2.0 * A_c * (Dr_beta_z / r) - 2.0 * third * sigma * A_lambda * div_beta
		+ FT2 * (2.0 * A_h * (DD_beta_rr / r) + 2.0 * A_lambda * Dr_beta_r)
		+ (1.0 - FT2) * (2.0 * A_a * (DD_beta_rr / r) + 2.0 * A_lambda * (beta_r / r));

	//     m     ^  i    ^  m       i              ^in    ^     m        ^  i       m    ^mn . .     i
	// beta  d Delta - Delta  d beta + (1/3)sigma (g   d (D beta ) + 2 Delta  D beta ) + g   D D beta
	//        m                m                        n  m                   m              m n
	*sDeltar += beta_r * DAr_Deltar + beta_z * DAz_Deltar - Deltar * Dr_beta_r - Deltaz * Dz_beta_r
		+ third * sigma * (g_a * Dr_div_beta + r * g_c * Dz_div_beta + 2.0 * Deltar * div_beta)
		+ g_a * Drr_beta_r + g_b * Dzz_beta_r + 2.0 * r * g_c * Drz_beta_r + g_h * DD_beta_rr;
	*sDeltaz += beta_r * DAr_Deltaz + beta_z * DAz_Deltaz - Deltar * Dr_beta_z - Deltaz * Dz_beta_z
		+ third * sigma * (r * g_c * Dr_div_beta + g_b * Dz_div_beta + 2.0 * Deltaz * div_beta)
		+ g_a * Drr_beta_z + g_b * Dzz_beta_z + 2.0 * r * g_c * Drz_beta_z + g_h * (Dr_beta_z / r);
}

void parabolic_shift(double *sbeta_r, double *sbeta_z, double *sdtbeta_r, double *sdtbeta_z,
	const double c, const double sDeltar, const double sDeltaz)
{
	*sbeta_r = c * sDeltar;
	*sbeta_z = c * sDeltaz;
	*sdtbeta_r = 0.0;
	*sdtbeta_z = 0.0;
}

void hyperbolic_shift(double *sbeta_r, double *sbeta_z, double *sdtbeta_r, double *sdtbeta_z,
	const double c, const double eta, const double sDeltar, const double sDeltaz,
	const double dtbeta_r, const double dtbeta_z)
{
	*sbeta_r = dtbeta_r;
	*sbeta_z = dtbeta_z;
	//*sbeta_r = c * dtbeta_r;
	//*sbeta_z = c * dtbeta_z;
	*sdtbeta_r = c * sDeltar - eta * dtbeta_r;
	*sdtbeta_z = c * sDeltaz - eta * dtbeta_z;
}

#include "tools.h"

#include "param.h"
#include "arrays.h"

void shift_sources(void)
{
	int k;

	if (strcmp(shift, "none") == 0)
	{
		#pragma omp parallel shared(sbeta_r, sbeta_z, sdtbeta_r, sdtbeta_z)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				// Zero-out components: (just to be sure).
				sbeta_r[k] = 0.0;
				sbeta_z[k] = 0.0;
				sdtbeta_r[k] = 0.0;
				sdtbeta_z[k] = 0.0;
			}
		}
	}
	else
	{
		#pragma omp parallel shared(sa, sb, sh, sc, slambda, sA_a, sA_b, sA_h, sA_c, sA_lambda, sphi, sK, sDeltar, sDeltaz)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				// Geometry modification.
				calculate_sources_geometry_shift(sa + k, sb + k, sh + k, sc + k, slambda + k,
					sA_a + k, sA_b + k, sA_h + k, sA_c + k, sA_lambda + k, sphi + k, sK + k,
					sDeltar + k, sDeltaz + k,
					r[k], a[k], b[k], h[k], c[k], lambda[k],
					A_a[k], A_b[k], A_h[k], A_c[k], A_lambda[k],
					beta_r[k], Dr_beta_r[k], Dz_beta_r[k], Drr_beta_r[k], Dzz_beta_r[k], Drz_beta_r[k],
					beta_z[k], Dr_beta_z[k], Dz_beta_z[k], Drr_beta_z[k], Dzz_beta_z[k], Drz_beta_z[k],
					DD_beta_rr[k], div_beta[k], Dr_div_beta[k], Dz_div_beta[k],
					DAr_a[k], DAr_b[k], DAr_h[k], DAr_c[k], DAr_lambda[k],
					DAz_a[k], DAz_b[k], DAz_h[k], DAz_c[k], DAz_lambda[k],
					DAr_A_a[k], DAr_A_b[k], DAr_A_h[k], DAr_A_c[k], DAr_A_lambda[k],
					DAz_A_a[k], DAz_A_b[k], DAz_A_h[k], DAz_A_c[k], DAz_A_lambda[k],
					DAr_phi[k], DAz_phi[k], DAr_K[k], DAz_K[k],
					g_a[k], g_b[k], g_h[k], g_c[k],
					Deltar[k], DAr_Deltar[k], DAz_Deltar[k],
					Deltaz[k], DAr_Deltaz[k], DAz_Deltaz[k],
					sigma);
			}
		}
		// Modification for lapse.
		if (!(strcmp(slicing, "maximal") == 0))
		{
			#pragma omp parallel shared(salpha)
			{
				#pragma omp for schedule(guided)
				for (k = 0; k < DIM; k++)
				{
					salpha[k] += beta_r[k] * DAr_alpha[k] + beta_z[k] * DAz_alpha[k];
				}
			}
		}

		// Type of shift.
		if (strcmp(shift, "parabolic") == 0)
		{
			#pragma omp parallel shared(sbeta_r, sbeta_z, sdtbeta_r, sdtbeta_z)
			{
				#pragma omp for schedule(guided)
				for (k = 0; k < DIM; k++)
				{
					parabolic_shift(sbeta_r + k, sbeta_z + k, sdtbeta_r + k, sdtbeta_z + k,
						driver_c, sDeltar[k], sDeltaz[k]);
				}
			}
		}
		else if (strcmp(shift, "hyperbolic") == 0)
		{
			#pragma omp parallel shared(sbeta_r, sbeta_z, sdtbeta_r, sdtbeta_z)
			{
				#pragma omp for schedule(guided)
				for (k = 0; k < DIM; k++)
				{
					hyperbolic_shift(sbeta_r + k, sbeta_z + k, sdtbeta_r + k, sdtbeta_z + k,
						driver_c, driver_eta, sDeltar[k], sDeltaz[k],
						dtbeta_r[k], dtbeta_z[k]);
				}
			}
		}
		else if (strcmp(shift, "static") == 0)
		{
			#pragma omp parallel shared(sbeta_r, sbeta_z, sdtbeta_r, sdtbeta_z)
			{
				#pragma omp for schedule(guided)
				for (k = 0; k < DIM; k++)
				{
					sbeta_r[k] = 0.;
					sbeta_z[k] = 0.;
					sdtbeta_r[k] = 0.;
					sdtbeta_z[k] = 0.;
				}
			}
		}
	}

	return;
}
