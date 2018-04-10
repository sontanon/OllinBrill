#include "tools.h"

#include "idata_brillwave.h"
#include "idata_schwarschild.h"
#include "idata_testbed.h"

#include "param.h"
#include "arrays.h"

void initial_lapse(void);
void initial_shift(void);

void initial(void)
{
	// Loop counters.
	int k;

	// Numbers.
	double aux, zero, one;

	zero = 0.0;
	one = 1.0;

	// Minkowski.
	// By default, initial date is always set to Minkowski first.
	#pragma omp parallel shared(alpha, K, a, b, h, c, Deltar, Deltaz, psi, phi, psi4)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			// Lapse.
			alpha[k] = one;

			// Extrinsic curvature.
			K[k] = zero;

			// Metric functions.
			a[k] = one;
			b[k] = one;
			h[k] = one;
			c[k] = zero;

			lambda[k] = zero;

			// Delta.
			Deltar[k] = zero;
			Deltaz[k] = zero;

			// Conformal factors.
			psi[k] = one;
			psi4[k] = one;
			phi[k] = zero;
		}
	}

	// Brill waves initial data.
	if (strcmp(idata, "BrillWave") == 0)
	{
		idata_brillwave();
	}
	// Schwarschild initial data.
	else if (strcmp(idata, "Schwarschild") == 0)
	{
		idata_schwarschild();
	}
	// Test-bed initial data.
	else if (strcmp(idata, "TestBed") == 0)
	{
		idata_testbed();
	}

	// Set initial lapse.
	initial_lapse();

#ifdef SHIFT
	// Set initial shift.
	if (!(strcmp(shift, "none") == 0))
	{
		initial_shift();
	}
#endif

	return;
}

// Set initial lapse.
void initial_lapse(void)
{
	int k;
	double aux;

	if (strcmp(ilapse, "one") == 0)
	{
		#pragma omp parallel shared(alpha)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				alpha[k] = 1.0;
			}
		}
	}
	else if (strcmp(ilapse, "isotropic") == 0)
	{
		#pragma omp parallel shared(alpha) private(aux)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				aux = 0.5 * schwar_mass / rr[k];
				alpha[k] = (1.0 - aux) / (1.0 + aux);
			}
		}
	}
	else if (strcmp(ilapse, "test") == 0)
	{
		#pragma omp parallel shared(alpha)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				alpha[k] = 1.0 + rr[k] * rr[k] * (exp(-rr[k] * rr[k])) / (1.0 + rr[k] * rr[k]);
			}
		}
	}
	else if (strcmp(ilapse, "precolapsed2") == 0)
	{
		#pragma omp parallel shared(alpha)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				alpha[k] = 1.0 / (psi[k] * psi[k]);
			}
		}
	}
	else if (strcmp(ilapse, "precolapsed4") == 0)
	{
		#pragma omp parallel shared(alpha)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				alpha[k] = 1.0 / (psi[k] * psi[k] * psi[k] * psi[k]);
			}
		}
	}
	else if (strcmp(ilapse, "gaussian") == 0)
	{
		#pragma omp parallel shared(alpha)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				alpha[k] = 1.0 + gl_alpha0 * rr[k] * rr[k] * (exp(-(rr[k] - gl_rr0) * (rr[k] - gl_rr0)) + exp(-(rr[k] + gl_rr0) * (rr[k] + gl_rr0))) / (1.0 + rr[k] * rr[k]);
			}
		}
	}

	return;
}

void initial_shift(void)
{
	int k;
	double aux;

	if (strcmp(ishift, "zero") == 0)
	{
		#pragma omp parallel shared(beta_r, beta_z)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				beta_r[k] = 0.0;
				beta_z[k] = 0.0;
			}
		}
	}
	else if (strcmp(ishift, "gaussian") == 0)
	{
		#pragma omp parallel shared(beta_r, beta_z)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				beta_r[k] = gs_beta0 * rr[k] * r[k] * (exp(-(rr[k] - gs_rr0) * (rr[k] - gs_rr0)) + exp(-(rr[k] + gs_rr0) * (rr[k] + gs_rr0))) / (1.0 + rr[k] * rr[k]);
				beta_z[k] = gs_beta0 * rr[k] * z[k] * (exp(-(rr[k] - gs_rr0) * (rr[k] - gs_rr0)) + exp(-(rr[k] + gs_rr0) * (rr[k] + gs_rr0))) / (1.0 + rr[k] * rr[k]);
			}
		}
	}

	return;
}

