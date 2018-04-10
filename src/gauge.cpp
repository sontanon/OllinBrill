#include "tools.h"

#include "param.h"
#include "arrays.h"
#include "derivatives.h"
#include "general_elliptic.h"
#include "maximal_slicing_elliptic_coefficients.h"

#include "geometry_alpha_derivatives.h"
#include "geometry_alpha_second_covariant_derivative.h"

#undef DEBUG
#undef DEBUG_ONE

// Forward declarations.
void maximal(void);
void one_plus_log(void);
void harmonic(void);
void isotropic(void);

void gauge(void)
{
	if (strcmp(slicing, "1+log") == 0)
	{
		// Call 1+log wrapper.
		one_plus_log();
	}
	else if (strcmp(slicing, "harmonic") == 0)
	{
		// Call harmonic wrapper.
		harmonic();
	}
	else if (strcmp(slicing, "maximal") == 0)
	{
		// Call maximal slicing wrapper.
		maximal();
	}
	else if (strcmp(slicing, "isotropic") == 0)
	{
		// Call isotropic slicing wrapper.
		isotropic();
	}

	return;
}

void isotropic(void)
{
	int k;

	// Now we calculate auxiliary-geometry quantities associated with the lapse.
	calculate_alpha_derivatives();
	calculate_alpha_second_covariant_derivative();

	return;
}

void one_plus_log(void)
{
	int k;

	// The gauge function is f(alpha) = gauge_f / alpha.
	#pragma omp parallel shared(falpha)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			falpha[k] = gauge_f / alpha[k];
		}
	}

	// Now we calculate auxiliary-geometry quantities associated with the lapse.
	calculate_alpha_derivatives();
	calculate_alpha_second_covariant_derivative();

	return;
}

void harmonic(void)
{
	int k;

	// The gauge function is f(alpha) = gauge_f.
	// Only necessary for initial time.
	#pragma omp parallel shared(falpha)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			falpha[k] = gauge_f;
		}
	}

	// Now we calculate auxiliary-geometry quantities associated with the lapse.
	calculate_alpha_derivatives();
	calculate_alpha_second_covariant_derivative();

	return;
}

void maximal(void)
{
	// To solve for maximal slicing we first have to construct the elliptic coefficients 
	// for the equation: D^2 alpha = alpha * (K^{ij} K_{ij}). 
	// In terms of the program arrays, the linear term is K2 and the Laplace operator will 
	// be written in terms of the derivatives of the physical metric.

	// Auxiliary variables.
	int k;

	// We only do an actual calculation for t > 0. At the initial time we have the lapse
	// tivially set to one and its derivatives to zero.
	if (t > 0.0)
	{
#ifdef VERBOSE
		printf("**************************************\n");
		printf("***                                ***\n");
		printf("***       Calculating maximal      ***\n");
		printf("***            slicing...          ***\n");
		printf("***                                ***\n");
#endif

		// Generate elliptic coefficients.
		#pragma omp parallel shared(ell_a, ell_b, ell_c, ell_d, ell_e, ell_f, ell_s)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				calculate_maximal_slicing_elliptic_coefficients(ell_a + k, ell_b + k, ell_c + k, ell_d + k, ell_e + k, ell_s + k, ell_f + k,
					r[k], a[k], b[k], h[k], c[k],
					Dr_a[k], Dr_b[k], Dr_h[k], Dr_c[k],
					Dz_a[k], Dz_b[k], Dz_h[k], Dz_c[k],
					Dr_phi[k], Dz_phi[k], psi4[k], A2[k]);
			}
		}

		// Call elliptic solver.
		general_elliptic(alpha, auxarray, ell_a, ell_b, ell_c, ell_d, ell_e, ell_s, ell_f, 1.0, nrobin, 1, 1, perm_flag, precond_flag, low_rank_flag);

		// Permutation.
		if (use_permutation)
		{
			// Increase preconditioner counter.
			perm_count += 1;

			// If permutation was just calculated, we can now use it.
			if (perm_flag == 2)
			{
				perm_flag = 1;
			}

			// Reset permutation flag every perm_count_max.
			if (perm_count % perm_count_max == 0)
			{
				// Setup to recalculate at next step.
				perm_flag = 2;
			}
		}
		// Preconditioner.
		// If we just did first solver, we can set preconditioner L flag.
		if (use_preconditioner && precond_flag == 0)
		{
			precond_flag = precond_L * 10 + 1;
		}

		// Low Rank Update.
		if (use_low_rank && low_rank_flag == 0)
		{
			low_rank_flag = 1;
		}

		// Now we calculate auxiliary-geometry quantities associated with the lapse.
		calculate_alpha_derivatives();
		calculate_alpha_second_covariant_derivative();

#ifdef DEBUG
		writeSingleFile(alpha, "ell_alpha.asc");
		writeSingleFile(Dr_alpha, "ell_Dr_alpha.asc");
		writeSingleFile(Dz_alpha, "ell_Dz_alpha.asc");
		writeSingleFile(Drr_alpha, "ell_Drr_alpha.asc");
		writeSingleFile(Dzz_alpha, "ell_Dzz_alpha.asc");
		writeSingleFile(Drz_alpha, "ell_Drz_alpha.asc");
		writeSingleFile(ell_a, "ell_a.asc");
		writeSingleFile(ell_b, "ell_b.asc");
		writeSingleFile(ell_c, "ell_c.asc");
		writeSingleFile(ell_d, "ell_d.asc");
		writeSingleFile(ell_e, "ell_e.asc");
		writeSingleFile(ell_s, "ell_s.asc");
		writeSingleFile(ell_f, "ell_f.asc");
#endif
	}
	// Set at t = 0 quantities.
	else if (t == 0.0)
	{
		#pragma omp parallel shared(alpha, falpha, Dr_alpha, Dz_alpha,\
		Drr_alpha, Drz_alpha, Dzz_alpha, \
		D2alpha_a, D2alpha_b, D2alpha_h, D2alpha_c, D2alpha_lambda, LaplaAlpha, DDalpha_r)
		{
		#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				alpha[k] = 1.0;
				falpha[k] = 1.0;
				Dr_alpha[k] = 0.0;
				Dz_alpha[k] = 0.0;
				Drr_alpha[k] = 0.0;
				Drz_alpha[k] = 0.0;
				Dzz_alpha[k] = 0.0;
				D2alpha_a[k] = 0.0;
				D2alpha_b[k] = 0.0;
				D2alpha_h[k] = 0.0;
				D2alpha_c[k] = 0.0;
				D2alpha_lambda[k] = 0.0;
				LaplaAlpha[k] = 0.0;
				DDalpha_r[k] = 0.0;
			}
		}
	}

	return;
}
