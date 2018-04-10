#include "tools.h"

#include "param.h"
#include "arrays.h"

void radbound(const double speed, const double varInf, const double *var, const double *Dr_var, const double *Dz_var, double *svar)
{
	MKL_INT i, j;

	double v_rr;

	// Boundaries on R.
	#pragma omp parallel shared(svar) private(v_rr)
	{
		#pragma omp for schedule(guided)
		for (j = 0; j < NzTotal - 1; j++)
		{
			v_rr = rr[IDX(NrTotal - 1, j)];

			svar[IDX(NrTotal - 1, j)] =
				-speed * (v_rr * Dr_var[IDX(NrTotal - 1, j)] / r[IDX(NrTotal - 1, j)]
					+ (var[IDX(NrTotal - 1, j)] - varInf) / v_rr);
		}
	}

	// Boundaries on Z.
	#pragma omp parallel shared(svar) private(v_rr)
	{
		#pragma omp for schedule(guided)
		for (i = 0; i < NrTotal - 1; i++)
		{
			v_rr = rr[IDX(i, NzTotal - 1)];

			svar[IDX(i, NzTotal - 1)] =
				-speed * (v_rr * Dz_var[IDX(i, NzTotal - 1)] / z[IDX(i, NzTotal - 1)]
					+ (var[IDX(i, NzTotal - 1)] - varInf) / v_rr);
		}
	}

	// Corner.
	double drr = sqrt(dr*dr + dz*dz);
	double Drr_var = (3.0 * var[IDX(NrTotal - 1, NzTotal - 1)] - 4.0 * var[IDX(NrTotal - 2, NzTotal - 2)] + var[IDX(NrTotal - 3, NzTotal - 3)]) / (2.0 * drr);

	v_rr = rr[IDX(NrTotal - 1, NzTotal - 1)];

	svar[IDX(NrTotal - 1, NzTotal - 1)] =
		-speed * (Drr_var + (var[IDX(NrTotal - 1, NzTotal - 1)] - varInf) / v_rr);
}
