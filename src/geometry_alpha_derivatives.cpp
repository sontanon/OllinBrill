#include "tools.h"

#include "param.h"
#include "arrays.h"

#include "derivatives.h"

void calculate_alpha_derivatives(void)
{
	int k;
	diff1r(Dr_alpha, alpha, 1);
	diff1z(Dz_alpha, alpha, 1);
	// Regularization.
#pragma omp parallel shared(auxarray)
	{
#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			auxarray[k] = Dr_alpha[k] / r[k];
		}
	}
	diff1r(DDalpha_r, auxarray, 1);

	diff2r(Drr_alpha, alpha, 1);
	diff2z(Dzz_alpha, alpha, 1);
	diff2rz(Drz_alpha, alpha, 1, 1);
}
