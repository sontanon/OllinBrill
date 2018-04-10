#include "tools.h"
#include "ah_param.h"

void calculate_constraints_full(double *ham, double *mom_r, double *mom_z,
	double *Clambda, double *CA_lambda,
	double *CDeltar, double *CDeltaz,
	const double r, const double a, const double h, const double lambda,
	const double A_a, const double A_h, const double A_lambda,
	const double g_a, const double g_b, const double g_c, const double g_lambda,
	const double Dr_g_a, const double Dz_g_b, const double Dr_g_c, const double Dz_g_c,
	const double R, const double A2, const double divA_r, const double divA_z,
	const double upA_a, const double upA_b, const double upA_c,
	const double Dr_phi, const double Dz_phi, const double psi4,
	const double Deltar, const double Deltaz,
	const double K, const double Dr_K, const double Dz_K,
	const double hdet, const double Dr_hdet, const double Dz_hdet)
{
	const double third = 1. / 3.;

	*ham = 0.5 * (R + 2. * third * K * K - A2);

	*mom_r = (divA_r + 6. * (upA_a * Dr_phi + r * upA_c * Dz_phi)
		- 2. * third * (g_a * Dr_K + r * g_c * Dz_K)) / psi4;
	*mom_z = (divA_z + 6. * (upA_b * Dz_phi + r * upA_c * Dr_phi)
		- 2. * third * (g_b * Dz_K + r * g_c * Dr_K)) / psi4;

	*CDeltar = Deltar - (-Dr_g_a - r * Dz_g_c - r * g_lambda
		- 0.5 * (g_a * Dr_hdet + r * g_c * Dz_hdet) / hdet);
	*CDeltaz = Deltaz - (-r * Dr_g_c - Dz_g_b - 2. * g_c
		- 0.5 * (r * g_c * Dr_hdet + g_b * Dz_hdet) / hdet);

	*Clambda = r * r * lambda - (a - h);
	*CA_lambda = r * r * A_lambda - (A_a - A_h);
}

void calculate_constraints_maximal(double *ham, double *mom_r, double *mom_z,
	double *Clambda, double *CA_lambda,
	double *CDeltar, double *CDeltaz,
	const double r, const double R, const double A2,
	const double divA_r, const double divA_z,
	const double upA_a, const double upA_b, const double upA_c,
	const double Dr_phi, const double Dz_phi, const double psi4,
	const double Deltar, const double Deltaz, const double Dr_g_a,
	const double Dz_g_c, const double Dr_g_c, const double Dz_g_b,
	const double g_a, const double g_b,
	const double g_c, const double g_lambda,
	const double lambda, const double a, const double h,
	const double A_lambda, const double A_a, const double A_h,
	const double hdet, const double Dr_hdet, const double Dz_hdet)
{
	const double third = 1. / 3.;

	// Simplify expressions if using maximal slicing.
	*ham = 0.5 * (R - A2);

	*mom_r = (divA_r + 6. * (upA_a * Dr_phi + r * upA_c * Dz_phi)) / psi4;
	*mom_z = (divA_z + 6. * (upA_b * Dz_phi + r * upA_c * Dr_phi)) / psi4;

	*CDeltar = Deltar - (-Dr_g_a - r * Dz_g_c - r * g_lambda
		- 0.5 * (g_a * Dr_hdet + r * g_c * Dz_hdet) / hdet);
	*CDeltaz = Deltaz - (-r * Dr_g_c - Dz_g_b - 2. * g_c
		- 0.5 * (r * g_c * Dr_hdet + g_b * Dz_hdet) / hdet);

	*Clambda = r * r * lambda - (a - h);
	*CA_lambda = r * r * A_lambda - (A_a - A_h);
}

#include "param.h"
#include "arrays.h"

#include "files.h"
#include "write_flags.h"

void constraints(const int l)
{
	int k;

	#pragma omp parallel shared(ham, mom_r, mom_z, \
	Clambda, CA_lambda, CDeltar, CDeltaz)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			calculate_constraints_full(ham + k, mom_r + k, mom_z + k,
				Clambda + k, CA_lambda + k,
				CDeltar + k, CDeltaz + k,
				r[k], a[k], h[k], lambda[k],
				A_a[k], A_h[k], A_lambda[k],
				g_a[k], g_b[k], g_c[k], g_lambda[k],
				Dr_g_a[k], Dz_g_b[k], Dr_g_c[k], Dz_g_c[k],
				RSCAL[k], A2[k], divA_r[k], divA_z[k],
				upA_a[k], upA_b[k], upA_c[k],
				Dr_phi[k], Dz_phi[k], psi4[k],
				Deltar[k], Deltaz[k],
				K[k], Dr_K[k], Dz_K[k],
				1.0, 0.0, 0.0);
//				hdet[k], Dr_hdet[k], Dz_hdet[k]);
		}
	}

	// FIND CONSTRAINTS NORM.
	double ham_norm = 0.0;
	double mom_r_norm = 0.0;
	double mom_z_norm = 0.0;
	double CDeltar_norm = 0.0;
	double CDeltaz_norm = 0.0;
	double Clambda_norm = 0.0;
	double CA_lambda_norm = 0.0;

	// Apply AH mask for t > ahfind_start_time.
#ifdef AH
	if (ahfind_flag && t >= ahfind_start_time && ahfind_mask_filter)
	{
#pragma omp parallel shared(ham, mom_r, mom_z, CDeltar, CDeltaz, Clambda, CA_lambda)
		{
#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				ham[k] *= (double)ah_mask_array[k];
				mom_r[k] *= (double)ah_mask_array[k];
				mom_z[k] *= (double)ah_mask_array[k];
				CDeltar[k] *= (double)ah_mask_array[k];
				CDeltaz[k] *= (double)ah_mask_array[k];
				Clambda[k] *= (double)ah_mask_array[k];
				CA_lambda[k] *= (double)ah_mask_array[k];
			}
		}
	}
#endif

	// Area element for norm calculation.
	double da = 1.0 / sqrt(DIM);

	// Loop integers.
	int i, j;

	// Calculate norms with boundary closing in.
	if (adaptive_boundary_norm)
	{
		if (l > NrTotal / 2)
		{
			printf("ERROR: Boundary advance has exceeded half of grid extension!\n");

			exit(-1);
		}


		#pragma omp parallel reduction(+:ham_norm, mom_r_norm, mom_z_norm, CDeltar_norm, CDeltaz_norm, Clambda_norm, CA_lambda_norm) private(j)
		{
			#pragma omp for schedule(guided)
			for (i = 0; i < NrTotal - l; i++)
			{
				for (j = 0; j < NzTotal - l; j++)
				{
					ham_norm += ham[IDX(i, j)] * ham[IDX(i, j)];
					mom_r_norm += mom_r[IDX(i, j)] * mom_r[IDX(i, j)];
					mom_z_norm += mom_z[IDX(i, j)] * mom_z[IDX(i, j)];
					CDeltar_norm += CDeltar[IDX(i, j)] * CDeltar[IDX(i, j)];
					CDeltaz_norm += CDeltaz[IDX(i, j)] * CDeltaz[IDX(i, j)];
					Clambda_norm += Clambda[IDX(i, j)] * Clambda[IDX(i, j)];
					CA_lambda_norm += CA_lambda[IDX(i, j)] * CA_lambda[IDX(i, j)];
				}
			}
		}

		da = 1.0 / sqrt((NrTotal - l) * (NzTotal - l));
		ham_norm = sqrt(ham_norm) * da;
		mom_r_norm = sqrt(mom_r_norm) * da;
		mom_z_norm = sqrt(mom_z_norm) * da;
		CDeltar_norm = sqrt(CDeltar_norm) * da;
		CDeltaz_norm = sqrt(CDeltaz_norm) * da;
		Clambda_norm = sqrt(Clambda_norm) * da;
		CA_lambda_norm = sqrt(CA_lambda_norm) * da;
	}
	else if (interior_norm)
	{
		#pragma omp parallel reduction(+:ham_norm, mom_r_norm, mom_z_norm, CDeltar_norm, CDeltaz_norm, Clambda_norm, CA_lambda_norm) private(j)
		{
			#pragma omp for schedule(guided)
			for (i = 0; i < NrTotal/2; i++)
			{
				for (j = 0; j < NzTotal/2; j++)
				{
					ham_norm += ham[IDX(i, j)] * ham[IDX(i, j)];
					mom_r_norm += mom_r[IDX(i, j)] * mom_r[IDX(i, j)];
					mom_z_norm += mom_z[IDX(i, j)] * mom_z[IDX(i, j)];
					CDeltar_norm += CDeltar[IDX(i, j)] * CDeltar[IDX(i, j)];
					CDeltaz_norm += CDeltaz[IDX(i, j)] * CDeltaz[IDX(i, j)];
					Clambda_norm += Clambda[IDX(i, j)] * Clambda[IDX(i, j)];
					CA_lambda_norm += CA_lambda[IDX(i, j)] * CA_lambda[IDX(i, j)];
				}
			}
		}

		da = 1.0 / sqrt((NrTotal/2) * (NzTotal/2));
		ham_norm = sqrt(ham_norm) * da;
		mom_r_norm = sqrt(mom_r_norm) * da;
		mom_z_norm = sqrt(mom_z_norm) * da;
		CDeltar_norm = sqrt(CDeltar_norm) * da;
		CDeltaz_norm = sqrt(CDeltaz_norm) * da;
		Clambda_norm = sqrt(Clambda_norm) * da;
		CA_lambda_norm = sqrt(CA_lambda_norm) * da;
	}
	else
	{
		ham_norm = cblas_dnrm2(DIM, ham, 1) * da;
		mom_r_norm = cblas_dnrm2(DIM, mom_r, 1) * da;
		mom_z_norm = cblas_dnrm2(DIM, mom_z, 1) * da;
		CDeltar_norm = cblas_dnrm2(DIM, CDeltar, 1) * da;
		CDeltaz_norm = cblas_dnrm2(DIM, CDeltar, 1) * da;
		Clambda_norm = cblas_dnrm2(DIM, Clambda, 1) * da;
		CA_lambda_norm = cblas_dnrm2(DIM, CA_lambda, 1) * da;
	}

	double tol;

	if (strcmp(order, "two") == 0)
	{
		tol = dr * dz;
	}
	else if (strcmp(order, "four") == 0)
	{
		tol = dr * dz * dr * dz;
		//tol = dr * dz;
	}

	// Print info to screen.
	printf("***                                ***\n");
	printf("***     Calculated constraints     ***\n");
	printf("***     in units of grid order:    ***\n");
	printf("***                                ***\n");
	printf("***     |ham|   =   %-6.3E      ***\n", ham_norm / tol);
	printf("***     |mom_r| =   %-6.3E      ***\n", mom_r_norm / tol);
	printf("***     |mom_z| =   %-6.3E      ***\n", mom_z_norm / tol);
	printf("***     |CDeltar|   = %-6.3E    ***\n", CDeltar_norm / tol);
	printf("***     |CDeltaz|   = %-6.3E    ***\n", CDeltaz_norm / tol);
	printf("***     |Clambda|   = %-6.3E    ***\n", Clambda_norm / tol);
	printf("***     |CA_lambda| = %-6.3E    ***\n", CA_lambda_norm / tol);
	printf("***                                ***\n");

	// Print to files.
	if (w_t_ham_norm)
		fprintf(f_t_ham_norm, "%9.18E\t%9.18E\n", t, ham_norm);
	if (w_t_mom_r_norm)
		fprintf(f_t_mom_r_norm, "%9.18E\t%9.18E\n", t, mom_r_norm);
	if (w_t_mom_z_norm)
		fprintf(f_t_mom_z_norm, "%9.18E\t%9.18E\n", t, mom_z_norm);
	if (w_t_CDeltar_norm)
		fprintf(f_t_CDeltar_norm, "%9.18E\t%9.18E\n", t, CDeltar_norm);
	if (w_t_CDeltaz_norm)
		fprintf(f_t_CDeltaz_norm, "%9.18E\t%9.18E\n", t, CDeltaz_norm);
	if (w_t_Clambda_norm)
		fprintf(f_t_Clambda_norm, "%9.18E\t%9.18E\n", t, Clambda_norm);
	if (w_t_CA_lambda_norm)
		fprintf(f_t_CA_lambda_norm, "%9.18E\t%9.18E\n", t, CA_lambda_norm);
}
