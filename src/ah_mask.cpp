#include "tools.h"

#include "param.h"
#include "arrays.h"

#include "ah_param.h"

double linear_interpolator(const int k1, const double k, const double *f)
{
	return f[k1] + (f[k1 + 1] - f[k1]) * (k - k1);
}

double cubic_interpolator(const int k1, const double k, const double *f, const double *df, const double dx)
{
	double f0 = f[k1];
	double f1 = f[k1 + 1];
	double df0 = dx * df[k1];
	double df1 = dx * df[k1 + 1];

	double a = 2.0 * f0 - 2.0 * f1 + df0 + df1;
	double b = -3.0 * f0 + 3.0 * f1 - 2.0 * df0 - df1;
	double c = df0;
	double d = f0;

	double x = k - (double)k1;

	return d + x * (c + x * (b + x * a));
}

void ah_mask(void)
{
	// Auxiliary integers.
	int i = 0, j = 0, k = 0;

	// Auxiliary doubles.
	double t_r, t_z, t_rr, t_th, t_k, t_R;

	// Set count to zero initially.
	ah_count = 0;

	// Check via cubic interpolation.
	if (ahfind_bicubic)
	{
		#pragma omp parallel shared(ah_mask_array) private(j, t_r, t_z, t_th, t_rr, t_k, k, t_R) 
		{
			#pragma omp for schedule(guided) \
			reduction(+:ah_count)
			for (i = ghost; i < NrTotal; i ++)
			{
				for (j = ghost; j < NzTotal; j++)
				{
					t_rr = rr[IDX(i, j)];

					// Check if we are trivially inside or outside the AH.
					if (t_rr > ah_max)
					{
						ah_mask_array[IDX(i, j)] = 1;
					}
					else if (t_rr < ah_min)
					{
						ah_mask_array[IDX(i, j)] = 0;
						// Increase count.
						ah_count += 1;
					}
					// Else we find values via interpolation.
					else
					{
						t_r = r[IDX(i, j)];
						t_z = z[IDX(i, j)];
						t_th = atan(t_r / t_z);
						// Find pair of angles between the obtained angle.
						t_k = t_th / dth;
						k = (int)floor(t_k);
						// Interpolate value cubicaly.
						t_R = cubic_interpolator(k, t_k, array_H, array_dH, dth);
						// Find if inside or outside AH.
						if (t_rr > t_R)
						{
							ah_mask_array[IDX(i ,j)] = 1;
						}
						else
						{
							ah_mask_array[IDX(i, j)] = 0;
							// Increase count.
							ah_count += 1;
						}
					}
				}
			}
		}

	}
	// Check via linear interpolation.
	else
	{
		#pragma omp parallel shared(ah_mask_array) private(j, t_r, t_z, t_th, t_rr, t_k, k, t_R) 
		{
			#pragma omp for schedule(guided) \
			reduction(+:ah_count)
			for (i = ghost; i < NrTotal; i ++)
			{
				for (j = ghost; j < NzTotal; j++)
				{
					t_rr = rr[IDX(i, j)];

					// Check if we are trivially inside or outside the AH.
					if (t_rr > ah_max)
					{
						ah_mask_array[IDX(i, j)] = 1;
					}
					else if (t_rr < ah_min)
					{
						ah_mask_array[IDX(i, j)] = 0;
						// Increase count.
						ah_count += 1;
					}
					// Else we find values via interpolation.
					else
					{
						t_r = r[IDX(i, j)];
						t_z = z[IDX(i, j)];
						t_th = atan(t_r / t_z);
						// Find pair of angles between the obtained angle.
						t_k = t_th / dth;
						k = (int)floor(t_k);
						// Interpolate value cubicaly.
						t_R = linear_interpolator(k, t_k, array_H);
						// Find if inside or outside AH.
						if (t_rr > t_R)
						{
							ah_mask_array[IDX(i ,j)] = 1;
						}
						else
						{
							ah_mask_array[IDX(i, j)] = 0;
							// Increase count.
							ah_count += 1;
						}
					}
				}
			}
		}
	}

	// Symmetries auxiliary count.
	int t_count;

	// Fill symmetries.
	for (k = 0; k < ghost; k++)
	{
		// Initialize t_count to zero,
		t_count = 0;

		// Loop in Z.
		#pragma omp parallel shared(ah_mask_array) 
		{
			#pragma omp for schedule(guided) \
			reduction(+:t_count)
			for (j = ghost - k; j < NzTotal; j++)
			{
				ah_mask_array[IDX(ghost - 1 - k, j)] = ah_mask_array[IDX(ghost + k, j)];
				// Increase count.
				if (ah_mask_array[IDX(ghost + k, j)] == 0)
				{
					t_count += 1;
				}
			}
		}
		
		// Update.
		ah_count += t_count;

		// Reset count.
		t_count = 0;

		// Loop in R.
		#pragma omp parallel shared(ah_mask_array)
		{
			#pragma omp for schedule(guided) \
			reduction(+:t_count)
			for (i = ghost - k; i < NrTotal; i++)
			{
				ah_mask_array[IDX(i, ghost - 1 - k)] = ah_mask_array[IDX(i, ghost + k)];
				// Increase count.
				if (ah_mask_array[IDX(i, ghost + k)] == 0)
				{
					t_count += 1;
				}
			}
		}

		// Update.
		ah_count += t_count;

		// Diagonal.
		ah_mask_array[IDX(ghost - 1 - k, ghost - 1 - k)] = ah_mask_array[IDX(ghost + k, ghost + k)];
		if (ah_mask_array[IDX(ghost + k, ghost + k)] == 0)
		{
			ah_count += 1;
		}
	}
}
