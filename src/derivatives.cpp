#include "param.h"
#include "tools.h"

void diff1r(double *dvar, const double *var, const int symr)
{
	// Auxiliary integers.
	int i, j;

	// Inverse spatial step.
	double idr = 1.0 / dr;

	// Constants.
	const double half = 0.5;
	const double twelfth = 1.0 / 12.0;

	// Second-order derivatives.
	if (strcmp(order, "two") == 0)
	{
#pragma omp parallel shared(dvar) private(j)
		{
			// Interior points.
#pragma omp for schedule(guided)
			for (i = ghost; i < NrTotal - 1; i++)
			{
				for (j = 0; j < NzTotal; j++)
				{
					dvar[IDX(i, j)] = idr * half * (var[IDX(i + 1, j)] - var[IDX(i - 1, j)]);
				}
			}
		}

		// Boundary.
#pragma omp parallel shared(dvar)
		{
#pragma omp for schedule(guided)
			for (j = 0; j < NzTotal; j++)
			{
				dvar[IDX(NrTotal - 1, j)] = idr * half * (3.0 * var[IDX(NrTotal - 1, j)] - 4.0 * var[IDX(NrTotal - 2, j)] + var[IDX(NrTotal - 3, j)]);
			}
		}
	}
	// Fourth-order derivatives.
	else if (strcmp(order, "four") == 0)
	{
#pragma omp parallel shared(dvar) private(j)
		{
			// Interior points.
#pragma omp for schedule(guided)
			for (i = ghost; i < NrTotal - 2; i++)
			{
				for (j = 0; j < NzTotal; j++)
				{
					dvar[IDX(i, j)] = idr * twelfth *  (8.0 * (var[IDX(i + 1, j)] - var[IDX(i - 1, j)]) - (var[IDX(i + 2, j)] - var[IDX(i - 2, j)]));
				}
			}
		}
#pragma omp parallel shared(dvar)
		{
			// Second-to-last and last points.
#pragma omp for schedule(guided)
			for (j = 0; j < NzTotal; j++)
			{
				dvar[IDX(NrTotal - 2, j)] = idr * twelfth * (3.0 * var[IDX(NrTotal - 1, j)] + 10.0 * var[IDX(NrTotal - 2, j)] - 18.0 * var[IDX(NrTotal - 3, j)] + 6.0 * var[IDX(NrTotal - 4, j)] - var[IDX(NrTotal - 5, j)]);
				dvar[IDX(NrTotal - 1, j)] = idr * twelfth * (25.0 * var[IDX(NrTotal - 1, j)] - 48.0 * var[IDX(NrTotal - 2, j)] + 36.0 * var[IDX(NrTotal - 3, j)] - 16.0 * var[IDX(NrTotal - 4, j)] + 3.0 * var[IDX(NrTotal - 5, j)]);
			}
		}
	}
	// Symmetries on axis.
	for (i = 0; i < ghost; i++)
	{
#pragma omp parallel shared(dvar)
		{
#pragma omp for schedule(guided)
			for (j = 0; j < NzTotal; j++)
			{
				dvar[IDX(ghost - 1 - i, j)] = -(double)(symr)* dvar[IDX(ghost + i, j)];
			}
		}
	}
}

void diff1z(double *dvar, const double *var, const int symz)
{
	int i, j;

	double idz = 1.0 / dz;

	const double half = 0.5;
	const double twelfth = 1.0 / 12.0;

	// Second-oder derivatives.
	if (strcmp(order, "two") == 0)
	{
		// Interior points.
#pragma omp parallel shared(dvar) private(i)
		{
#pragma omp for schedule(guided)
			for (j = ghost; j < NzTotal - 1; j++)
			{
				for (i = 0; i < NrTotal; i++)
				{
					dvar[IDX(i, j)] = idz * half * (var[IDX(i, j + 1)] - var[IDX(i, j - 1)]);
				}
			}
		}

		// Boundary.
#pragma omp parallel shared(dvar)
		{
#pragma omp for schedule(guided)
			for (i = 0; i < NrTotal; i++)
			{
				dvar[IDX(i, NzTotal - 1)] = idz * half * (3.0 * var[IDX(i, NzTotal - 1)] - 4.0 * var[IDX(i, NzTotal - 2)] + var[IDX(i, NzTotal - 3)]);
			}
		}
	}
	// Fourth-order derivatives.
	else if (strcmp(order, "four") == 0)
	{
#pragma omp parallel shared(dvar) private(i)
		{
			// Interior points.
#pragma omp for schedule(guided)
			for (j = ghost; j < NzTotal - 2; j++)
			{
				for (i = 0; i < NrTotal; i++)
				{
					dvar[IDX(i, j)] = idz * twelfth * (8.0 * (var[IDX(i, j + 1)] - var[IDX(i, j - 1)]) - (var[IDX(i, j + 2)] - var[IDX(i, j - 2)]));
				}
			}
		}

#pragma omp parallel shared(dvar)
		{
			// Second-to-last and last points.
#pragma omp for schedule(guided)
			for (i = 0; i < NrTotal; i++)
			{
				dvar[IDX(i, NzTotal - 2)] = idz * twelfth * (3.0 * var[IDX(i, NzTotal - 1)] + 10.0 * var[IDX(i, NzTotal - 2)] - 18.0 * var[IDX(i, NzTotal - 3)] + 6.0 * var[IDX(i, NzTotal - 4)] - var[IDX(i, NzTotal - 5)]);
				dvar[IDX(i, NzTotal - 1)] = idz * twelfth * (25.0 * var[IDX(i, NzTotal - 1)] - 48.0 * var[IDX(i, NzTotal - 2)] + 36.0 * var[IDX(i, NzTotal - 3)] - 16.0 * var[IDX(i, NzTotal - 4)] + 3.0 * var[IDX(i, NzTotal - 5)]);
			}
		}
	}

	// Symmetries on equator.
	for (j = 0; j < ghost; j++)
	{
#pragma omp parallel shared(dvar)
		{
#pragma omp for schedule(guided)
			for (i = 0; i < NrTotal; i++)
			{
				dvar[IDX(i, ghost - 1 - j)] = -(double)(symz)* dvar[IDX(i, ghost + j)];
			}
		}
	}
}

void diff2r(double *dvar, const double *var, const int symr)
{
	int i, j;

	double idr2 = 1.0 / (dr * dr);

	const double twelfth = 1.0 / 12.0;

	// Second-order derivatives.
	if (strcmp(order, "two") == 0)
	{
#pragma omp parallel shared(dvar) private(j)
		{
			// Interior points.
#pragma omp for schedule(guided)
			for (i = ghost; i < NrTotal - 1; i++)
			{
				for (j = 0; j < NzTotal; j++)
				{
					dvar[IDX(i, j)] = idr2 * (var[IDX(i + 1, j)] - 2.0 * var[IDX(i, j)] + var[IDX(i - 1, j)]);
				}
			}
		}

		// Boundary points.
#pragma omp parallel shared(dvar)
		{
#pragma omp for schedule(guided)
			for (j = 0; j < NzTotal; j++)
			{
				dvar[IDX(NrTotal - 1, j)] = idr2 * (2.0 * var[IDX(NrTotal - 1, j)] - 5.0 * var[IDX(NrTotal - 2, j)] + 4.0 * var[IDX(NrTotal - 3, j)] - var[IDX(NrTotal - 4, j)]);

			}
		}

	}
	// Fourth-order derivatives.
	else if (strcmp(order, "four") == 0)
	{
#pragma omp parallel shared(dvar) private(j)
		{
			// Interior points.
#pragma omp for schedule(guided)
			for (i = ghost; i < NrTotal - 2; i++)
			{
				for (j = 0; j < NzTotal; j++)
				{
					dvar[IDX(i, j)] = -idr2 * twelfth * (30.0 * var[IDX(i, j)] - 16.0 * (var[IDX(i + 1, j)] + var[IDX(i - 1, j)]) + var[IDX(i + 2, j)] + var[IDX(i - 2, j)]);
				}
			}
		}

		// Second-to-last and last point.
#pragma omp parallel shared(dvar)
		{
#pragma omp for schedule(guided)
			for (j = 0; j < NzTotal; j++)
			{
				dvar[IDX(NrTotal - 2, j)] = idr2 * twelfth * (10.0 * var[IDX(NrTotal - 1, j)] - 15.0 * var[IDX(NrTotal - 2, j)] - 4.0 * var[IDX(NrTotal - 3, j)] + 14.0 * var[IDX(NrTotal - 4, j)] - 6.0 * var[IDX(NrTotal - 5, j)] + var[IDX(NrTotal - 6, j)]);

				dvar[IDX(NrTotal - 1, j)] = idr2 * twelfth * (45.0 * var[IDX(NrTotal - 1, j)] - 154.0 * var[IDX(NrTotal - 2, j)] + 214.0 * var[IDX(NrTotal - 3, j)] - 156.0 * var[IDX(NrTotal - 4, j)] + 61.0 * var[IDX(NrTotal - 5, j)] - 10.0 * var[IDX(NrTotal - 6, j)]);
			}
		}
	}

	// Symmetries on axis.
	for (i = 0; i < ghost; i++)
	{
#pragma omp parallel shared(dvar)
		{
#pragma omp for schedule(guided)
			for (j = 0; j < NzTotal; j++)
			{
				dvar[IDX(ghost - 1 - i, j)] = (double)(symr)* dvar[IDX(ghost + i, j)];
			}
		}
	}
}

void diff2z(double *dvar, const double *var, const int symz)
{
	int i, j;

	double idz2 = 1.0 / (dz * dz);

	const double twelfth = 1.0 / 12.0;

	// Second-order derivatives.
	if (strcmp(order, "two") == 0)
	{
#pragma omp parallel shared(dvar) private(i)
		{
			// Interior points.
#pragma omp for schedule(guided)
			for (j = ghost; j < NzTotal - 1; j++)
			{
				for (i = 0; i < NrTotal; i++)
				{
					dvar[IDX(i, j)] = idz2 * (var[IDX(i, j + 1)] - 2.0 * var[IDX(i, j)] + var[IDX(i, j - 1)]);
				}
			}
		}

		// Boundary points.
#pragma omp parallel shared(dvar)
		{
#pragma omp for schedule(guided)
			for (i = 0; i < NrTotal; i++)
			{
				dvar[IDX(i, NzTotal - 1)] = idz2 * (2.0 * var[IDX(i, NzTotal - 1)] - 5.0 * var[IDX(i, NzTotal - 2)] + 4.0 * var[IDX(i, NzTotal - 3)] - var[IDX(i, NzTotal - 4)]);
			}
		}

	}
	// Fourth-order derivatives.
	else if (strcmp(order, "four") == 0)
	{
#pragma omp parallel shared(dvar) private(i)
		{
			// Interior points.
#pragma omp for schedule(guided)
			for (j = ghost; j < NzTotal - 2; j++)
			{
				for (i = 0; i < NrTotal; i++)
				{
					dvar[IDX(i, j)] = -idz2 * twelfth * (30.0 * var[IDX(i, j)] - 16.0 * (var[IDX(i, j + 1)] + var[IDX(i, j - 1)]) + var[IDX(i, j + 2)] + var[IDX(i, j - 2)]);
				}
			}
		}

		// Second-to-last and last point.
#pragma omp parallel shared(dvar)
		{
#pragma omp for schedule(guided)
			for (i = 0; i < NrTotal; i++)
			{
				dvar[IDX(i, NzTotal - 2)] = idz2 * twelfth * (10.0 * var[IDX(i, NzTotal - 1)] - 15.0 * var[IDX(i, NzTotal - 2)] - 4.0 * var[IDX(i, NzTotal - 3)] + 14.0 * var[IDX(i, NzTotal - 4)] - 6.0 * var[IDX(i, NzTotal - 5)] + var[IDX(i, NzTotal - 6)]);

				dvar[IDX(i, NzTotal - 1)] = idz2 * twelfth * (45.0 * var[IDX(i, NzTotal - 1)] - 154.0 * var[IDX(i, NzTotal - 2)] + 214.0 * var[IDX(i, NzTotal - 3)] - 156.0 * var[IDX(i, NzTotal - 4)] + 61.0 * var[IDX(i, NzTotal - 5)] - 10.0 * var[IDX(i, NzTotal - 6)]);
			}
		}
	}

	// Symmetries on equator.
	for (j = 0; j < ghost; j++)
	{
#pragma omp parallel shared(dvar)
		{
#pragma omp for schedule(guided)
			for (i = 0; i < NrTotal; i++)
			{
				dvar[IDX(i, ghost - 1 - j)] = (double)(symz)*dvar[IDX(i, ghost + j)];
			}
		}
	}
}

void diff2rz(double *dvar, const double *var, const int symr, const int symz)
{
	int i, j;

	double idr = 1.0 / dr, idz = 1.0 / dz;
	double idrz = idr * idz;

	const double quarter = 0.25;
	const double i48 = 1.0 / 48.0;
	const double i144 = 1.0 / 144.0;

	// Second-order derivatives.
	if (strcmp(order, "two") == 0)
	{
		// Interior points.
#pragma omp parallel shared(dvar) private(j)
		{
#pragma omp for schedule(guided)
			for (i = ghost; i < NrTotal - 1; i++)
			{
				for (j = ghost; j < NzTotal - 1; j++)
				{
					dvar[IDX(i, j)] = idrz * quarter * (var[IDX(i + 1, j + 1)] + var[IDX(i - 1, j - 1)] - var[IDX(i + 1, j - 1)] - var[IDX(i - 1, j + 1)]);
				}

				// Last point on z.
				dvar[IDX(i, NzTotal - 1)] = idrz * quarter * ((3.0 * var[IDX(i + 1, NzTotal - 1)] - 4.0 * var[IDX(i + 1, NzTotal - 2)] + var[IDX(i + 1, NzTotal - 3)]) - (3.0 * var[IDX(i - 1, NzTotal - 1)] - 4.0 * var[IDX(i - 1, NzTotal - 2)] + var[IDX(i - 1, NzTotal - 3)]));
			}
		}

		// Last point on r
#pragma omp parallel shared(dvar)
		{
#pragma omp for schedule(guided)
			for (j = ghost; j < NzTotal - 1; j++)
			{
				dvar[IDX(NrTotal - 1, j)] = idrz * quarter * ((3.0 * var[IDX(NrTotal - 1, j + 1)] - 4.0 * var[IDX(NrTotal - 2, j + 1)] + var[IDX(NrTotal - 3, j + 1)]) - (3.0 * var[IDX(NrTotal - 1, j - 1)] - 4.0 * var[IDX(NrTotal - 2, j - 1)] + var[IDX(NrTotal - 3, j - 1)]));
			}
		}

		// Corner.
		dvar[IDX(NrTotal - 1, NzTotal - 1)] = idrz * quarter * (9.0 * var[IDX(NrTotal - 1, NzTotal - 1)] + 16.0 * var[IDX(NrTotal - 2, NzTotal - 2)] + var[IDX(NrTotal - 3, NzTotal - 3)] - 12.0 * (var[IDX(NrTotal - 2, NzTotal - 1)] + var[IDX(NrTotal - 1, NzTotal - 2)]) + 3.0 * (var[IDX(NrTotal - 1, NzTotal - 3)] + var[IDX(NrTotal - 3, NzTotal - 1)]) - 4.0 * (var[IDX(NrTotal - 2, NzTotal - 3)] + var[IDX(NrTotal - 3, NzTotal - 2)]));
	}
	// Fourth-order derivatives.
	else if (strcmp(order, "four") == 0)
	{
		// Interior points.
#pragma omp parallel shared(dvar) private(j)
		{
#pragma omp for schedule(guided)
			for (i = ghost; i < NrTotal - 2; i++)
			{
				for (j = ghost; j < NzTotal - 2; j++)
				{
					dvar[IDX(i, j)] = idrz * i48 * (16.0 * (var[IDX(i + 1, j + 1)] + var[IDX(i - 1, j - 1)] - var[IDX(i + 1, j - 1)] - var[IDX(i - 1, j + 1)]) - (var[IDX(i + 2, j + 2)] + var[IDX(i - 2, j - 2)] - var[IDX(i + 2, j - 2)] - var[IDX(i - 2, j + 2)]));
				}

				// Second-to-last and last points in Z.
				dvar[IDX(i, NzTotal - 2)] = -idrz * i144 * ((var[IDX(i - 2, NzTotal - 5)] - var[IDX(i + 2, NzTotal - 5)]) - 6.0 * (var[IDX(i - 2, NzTotal - 4)] - var[IDX(i + 2, NzTotal - 4)]) + 18.0 * (var[IDX(i - 2, NzTotal - 3)] - var[IDX(i + 2, NzTotal - 3)]) - 10.0 * (var[IDX(i - 2, NzTotal - 2)] - var[IDX(i + 2, NzTotal - 2)]) - 3.0 * (var[IDX(i - 2, NzTotal - 1)] - var[IDX(i + 2, NzTotal - 1)]) - 8.0 * (var[IDX(i - 1, NzTotal - 5)] - var[IDX(i + 1, NzTotal - 5)]) + 48.0 * (var[IDX(i - 1, NzTotal - 4)] - var[IDX(i + 1, NzTotal - 4)]) - 144.0 * (var[IDX(i - 1, NzTotal - 3)] - var[IDX(i + 1, NzTotal - 3)]) + 80.0 * (var[IDX(i - 1, NzTotal - 2)] - var[IDX(i + 1, NzTotal - 2)]) + 24.0 * (var[IDX(i - 1, NzTotal - 1)] - var[IDX(i + 1, NzTotal - 1)]));

				dvar[IDX(i, NzTotal - 1)] = -idrz * i144 * (-3.0 * (var[IDX(i - 2, NzTotal - 5)] - var[IDX(i + 2, NzTotal - 5)]) + 16.0 * (var[IDX(i - 2, NzTotal - 4)] - var[IDX(i + 2, NzTotal - 4)]) - 36.0 * (var[IDX(i - 2, NzTotal - 3)] - var[IDX(i + 2, NzTotal - 3)]) + 48.0 * (var[IDX(i - 2, NzTotal - 2)] - var[IDX(i + 2, NzTotal - 2)]) - 25.0 * (var[IDX(i - 2, NzTotal - 1)] - var[IDX(i + 2, NzTotal - 1)]) + 24.0 * (var[IDX(i - 1, NzTotal - 5)] - var[IDX(i + 1, NzTotal - 5)]) - 128.0 * (var[IDX(i - 1, NzTotal - 4)] - var[IDX(i + 1, NzTotal - 4)]) + 288.0 * (var[IDX(i - 1, NzTotal - 3)] - var[IDX(i + 1, NzTotal - 3)]) - 384.0 * (var[IDX(i - 1, NzTotal - 2)] - var[IDX(i + 1, NzTotal - 2)]) + 200.0 * (var[IDX(i - 1, NzTotal - 1)] - var[IDX(i + 1, NzTotal - 1)]));


			}

		}

		// Second-to-last and last points in R.
#pragma omp parallel shared(dvar)
		{
#pragma omp for schedule(guided)
			for (j = ghost; j < NzTotal - 2; j++)
			{
				dvar[IDX(NrTotal - 2, j)] = -idrz * i144 * ((var[IDX(NrTotal - 5, j - 2)] - var[IDX(NrTotal - 5, j + 2)]) - 6.0 * (var[IDX(NrTotal - 4, j - 2)] - var[IDX(NrTotal - 4, j + 2)]) + 18.0 * (var[IDX(NrTotal - 3, j - 2)] - var[IDX(NrTotal - 3, j + 2)]) - 10.0 * (var[IDX(NrTotal - 2, j - 2)] - var[IDX(NrTotal - 2, j + 2)]) - 3.0 * (var[IDX(NrTotal - 1, j - 2)] - var[IDX(NrTotal - 1, j + 2)]) - 8.0 * (var[IDX(NrTotal - 5, j - 1)] - var[IDX(NrTotal - 5, j + 1)]) + 48.0 * (var[IDX(NrTotal - 4, j - 1)] - var[IDX(NrTotal - 4, j + 1)]) - 144.0 * (var[IDX(NrTotal - 3, j - 1)] - var[IDX(NrTotal - 3, j + 1)]) + 80.0 * (var[IDX(NrTotal - 2, j - 1)] - var[IDX(NrTotal - 2, j + 1)]) + 24.0 * (var[IDX(NrTotal - 1, j - 1)] - var[IDX(NrTotal - 1, j + 1)]));

				dvar[IDX(NrTotal - 1, j)] = -idrz * i144 * (-3.0 * (var[IDX(NrTotal - 5, j - 2)] - var[IDX(NrTotal - 5, j + 2)]) + 16.0 * (var[IDX(NrTotal - 4, j - 2)] - var[IDX(NrTotal - 4, j + 2)]) - 36.0 * (var[IDX(NrTotal - 3, j - 2)] - var[IDX(NrTotal - 3, j + 2)]) + 48.0 * (var[IDX(NrTotal - 2, j - 2)] - var[IDX(NrTotal - 2, j + 2)]) - 25.0 * (var[IDX(NrTotal - 1, j - 2)] - var[IDX(NrTotal - 1, j + 2)]) + 24.0 * (var[IDX(NrTotal - 5, j - 1)] - var[IDX(NrTotal - 5, j + 1)]) - 128.0 * (var[IDX(NrTotal - 4, j - 1)] - var[IDX(NrTotal - 4, j + 1)]) + 288.0 * (var[IDX(NrTotal - 3, j - 1)] - var[IDX(NrTotal - 3, j + 1)]) - 384.0 * (var[IDX(NrTotal - 2, j - 1)] - var[IDX(NrTotal - 2, j + 1)]) + 200.0 * (var[IDX(NrTotal - 1, j - 1)] - var[IDX(NrTotal - 1, j + 1)]));
			}
		}

		// Corner: 4 points.
		dvar[IDX(NrTotal - 2, NzTotal - 2)] = idrz * i144 * (var[IDX(NrTotal - 5, NzTotal - 5)] + 36.0 * var[IDX(NrTotal - 4, NzTotal - 4)] + 324.0 * var[IDX(NrTotal - 3, NzTotal - 3)] + 100.0 * var[IDX(NrTotal - 2, NzTotal - 2)] + 9 * var[IDX(NrTotal - 1, NzTotal - 1)] - 6.0 * (var[IDX(NrTotal - 5, NzTotal - 4)] + var[IDX(NrTotal - 4, NzTotal - 5)]) + 18.0 * (var[IDX(NrTotal - 5, NzTotal - 3)] + var[IDX(NrTotal - 3, NzTotal - 5)]) - 10.0 * (var[IDX(NrTotal - 5, NzTotal - 2)] + var[IDX(NrTotal - 2, NzTotal - 5)]) - 3.0 * (var[IDX(NrTotal - 5, NzTotal - 1)] + var[IDX(NrTotal - 1, NzTotal - 5)]) - 108.0 * (var[IDX(NrTotal - 4, NzTotal - 3)] + var[IDX(NrTotal - 3, NzTotal - 4)]) + 60.0 * (var[IDX(NrTotal - 4, NzTotal - 2)] + var[IDX(NrTotal - 2, NzTotal - 4)]) + 18.0 * (var[IDX(NrTotal - 4, NzTotal - 1)] + var[IDX(NrTotal - 1, NzTotal - 4)]) - 180.0 * (var[IDX(NrTotal - 3, NzTotal - 2)] + var[IDX(NrTotal - 2, NzTotal - 3)]) - 54.0 * (var[IDX(NrTotal - 3, NzTotal - 1)] + var[IDX(NrTotal - 1, NzTotal - 3)]) + 30.0 * (var[IDX(NrTotal - 2, NzTotal - 1)] + var[IDX(NrTotal - 1, NzTotal - 2)]));

		dvar[IDX(NrTotal - 2, NzTotal - 1)] = idrz * i144 * (-3.0 * var[IDX(NrTotal - 5, NzTotal - 5)] - 96.0 * var[IDX(NrTotal - 4, NzTotal - 4)] - 648.0 * var[IDX(NrTotal - 3, NzTotal - 3)] - 480.0 * var[IDX(NrTotal - 2, NzTotal - 2)] + 75.0 * var[IDX(NrTotal - 1, NzTotal - 1)] + 16.0 * var[IDX(NrTotal - 5, NzTotal - 4)] + 18.0 * var[IDX(NrTotal - 4, NzTotal - 5)] - 36.0 * var[IDX(NrTotal - 5, NzTotal - 3)] - 54.0 * var[IDX(NrTotal - 3, NzTotal - 5)] + 48.0 * var[IDX(NrTotal - 5, NzTotal - 2)] + 30.0 * var[IDX(NrTotal - 2, NzTotal - 5)] - 25.0 * var[IDX(NrTotal - 5, NzTotal - 1)] + 9.0 * var[IDX(NrTotal - 1, NzTotal - 5)] + 216.0 * var[IDX(NrTotal - 4, NzTotal - 3)] + 288.0 * var[IDX(NrTotal - 3, NzTotal - 4)] - 288.0 * var[IDX(NrTotal - 4, NzTotal - 2)] - 160.0 * var[IDX(NrTotal - 2, NzTotal - 4)] + 150.0 * var[IDX(NrTotal - 4, NzTotal - 1)] - 48.0 * var[IDX(NrTotal - 1, NzTotal - 4)] + 864.0 * var[IDX(NrTotal - 3, NzTotal - 2)] + 360.0 * var[IDX(NrTotal - 2, NzTotal - 3)] - 450.0 * var[IDX(NrTotal - 3, NzTotal - 1)] + 108.0 * var[IDX(NrTotal - 1, NzTotal - 3)] + 250.0 * var[IDX(NrTotal - 2, NzTotal - 1)] - 144.0 * var[IDX(NrTotal - 1, NzTotal - 2)]);

		dvar[IDX(NrTotal - 1, NzTotal - 2)] = idrz * i144 * (-3.0 * var[IDX(NrTotal - 5, NzTotal - 5)] - 96.0 * var[IDX(NrTotal - 4, NzTotal - 4)] - 648.0 * var[IDX(NrTotal - 3, NzTotal - 3)] - 480.0 * var[IDX(NrTotal - 2, NzTotal - 2)] + 75.0 * var[IDX(NrTotal - 1, NzTotal - 1)] + 18.0 * var[IDX(NrTotal - 5, NzTotal - 4)] + 16.0 * var[IDX(NrTotal - 4, NzTotal - 5)] - 54.0 * var[IDX(NrTotal - 5, NzTotal - 3)] - 36.0 * var[IDX(NrTotal - 3, NzTotal - 5)] + 30.0 * var[IDX(NrTotal - 5, NzTotal - 2)] + 48.0 * var[IDX(NrTotal - 2, NzTotal - 5)] + 9.0 * var[IDX(NrTotal - 5, NzTotal - 1)] - 25.0 * var[IDX(NrTotal - 1, NzTotal - 5)] + 288.0 * var[IDX(NrTotal - 4, NzTotal - 3)] + 216.0 * var[IDX(NrTotal - 3, NzTotal - 4)] - 160.0 * var[IDX(NrTotal - 4, NzTotal - 2)] - 288.0 * var[IDX(NrTotal - 2, NzTotal - 4)] - 48.0 * var[IDX(NrTotal - 4, NzTotal - 1)] + 150.0 * var[IDX(NrTotal - 1, NzTotal - 4)] + 360.0 * var[IDX(NrTotal - 3, NzTotal - 2)] + 864.0 * var[IDX(NrTotal - 2, NzTotal - 3)] + 108.0 * var[IDX(NrTotal - 3, NzTotal - 1)] - 450.0 * var[IDX(NrTotal - 1, NzTotal - 3)] - 144.0 * var[IDX(NrTotal - 2, NzTotal - 1)] + 250.0 * var[IDX(NrTotal - 1, NzTotal - 2)]);

		dvar[IDX(NrTotal - 1, NzTotal - 1)] = idrz * i144 * (9.0 * var[IDX(NrTotal - 5, NzTotal - 5)] + 256.0 * var[IDX(NrTotal - 4, NzTotal - 4)] + 1296.0 * var[IDX(NrTotal - 3, NzTotal - 3)] + 2304.0 * var[IDX(NrTotal - 2, NzTotal - 2)] + 625.0 * var[IDX(NrTotal - 1, NzTotal - 1)] - 48.0 * (var[IDX(NrTotal - 5, NzTotal - 4)] + var[IDX(NrTotal - 4, NzTotal - 5)]) + 108.0 * (var[IDX(NrTotal - 5, NzTotal - 3)] + var[IDX(NrTotal - 3, NzTotal - 5)]) - 144.0 * (var[IDX(NrTotal - 5, NzTotal - 2)] + var[IDX(NrTotal - 2, NzTotal - 5)]) + 75.0 * (var[IDX(NrTotal - 5, NzTotal - 1)] + var[IDX(NrTotal - 1, NzTotal - 5)]) - 576.0 * (var[IDX(NrTotal - 4, NzTotal - 3)] + var[IDX(NrTotal - 3, NzTotal - 4)]) + 768.0 * (var[IDX(NrTotal - 4, NzTotal - 2)] + var[IDX(NrTotal - 2, NzTotal - 4)]) - 400.0 * (var[IDX(NrTotal - 4, NzTotal - 1)] + var[IDX(NrTotal - 1, NzTotal - 4)]) - 1728.0 * (var[IDX(NrTotal - 3, NzTotal - 2)] + var[IDX(NrTotal - 2, NzTotal - 3)]) + 900.0 * (var[IDX(NrTotal - 3, NzTotal - 1)] + var[IDX(NrTotal - 1, NzTotal - 3)]) - 1200.0 * (var[IDX(NrTotal - 2, NzTotal - 1)] + var[IDX(NrTotal - 1, NzTotal - 2)]));

	}

	// Symmetries on axis and equator.
	int k;
	for (k = 0; k < ghost; k++)
	{
#pragma omp parallel shared(dvar)
		{
#pragma omp for schedule(guided)
			for (i = ghost - k; i < NrTotal; i++)
			{
				dvar[IDX(i, ghost - 1 - k)] = -(double)(symz)*dvar[IDX(i, ghost + k)];
			}
		}

#pragma omp parallel shared (dvar)
		{
#pragma omp for schedule(guided)
			for (j = ghost - k; j < NzTotal; j++)
			{
				dvar[IDX(ghost - 1 - k, j)] = -(double)(symr)*dvar[IDX(ghost + k, j)];
			}
		}
		// Corner.
		dvar[IDX(ghost - 1 - k, ghost - 1 - k)] = (double)(symr*symz)*dvar[IDX(ghost + k, ghost + k)];
	}
}
