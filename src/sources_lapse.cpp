#include "tools.h"

#include "param.h"
#include "arrays.h"

void calculate_sources_lapse(double *salpha, const double falpha, const double alpha, const double K)
{
	// Bona-Masso family.
	*salpha = -falpha * alpha * alpha * K;

	return;
}

void sources_lapse(void)
{
	int k;

	if (strcmp(slicing, "maximal") == 0)
	{
		#pragma omp parallel shared(salpha)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				salpha[k] = 0.0;
			}
		}
	}
	else
	{
		#pragma omp parallel shared(salpha)
		{
			#pragma omp for schedule(guided)
			for (k = 0; k < DIM; k++)
			{
				calculate_sources_lapse(salpha + k, falpha[k], alpha[k], K[k]);
			}
		}
	}

	return;
}
