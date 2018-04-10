#include "tools.h"

#include "param.h"

#define MAIN_FILE
#include "ah_param.h"

void ah_init(void)
{
	printf("**************************************\n");
	printf("***                                ***\n");
	printf("***     Initializing AH finder...  ***\n");
	printf("***                                ***\n");

	// The first thing we need to do is to initialize arrays in which we will store
	// the values of H and dH. We also need arrays for a linspace of theta.
	// Our arrays will have a NthTotal = min(NrTotal, NzTotal) points.
	NthTotal = MIN(NrInterior, NzInterior) + 1;

	// Allocate memory.
	array_th = (double *)malloc(NthTotal * sizeof(double));
	array_H  = (double *)malloc(NthTotal * sizeof(double));
	array_dH = (double *)malloc(NthTotal * sizeof(double));
	array_a = (double *)malloc(NthTotal * sizeof(double));
	array_b = (double *)malloc(NthTotal * sizeof(double));
	array_h = (double *)malloc(NthTotal * sizeof(double));
	array_c = (double *)malloc(NthTotal * sizeof(double));
	array_phi = (double *)malloc(NthTotal * sizeof(double));

	// Get step sizes.
	drr = 0.25 * sqrt(dr * dz);
	dth = M_PI / (2. * (NthTotal - 1));

	// Fill theta linspace.
	int k = 0;

	for (k = 0; k < NthTotal; k++)
	{
		array_th[k] = k * dth;
	}

	// AH mask.
	ah_mask_array = (int *)malloc(DIM * sizeof(int));
	// Set mask to 1.
	#pragma omp parallel shared(ah_mask_array)
	{
		#pragma omp for schedule(guided) 
		for (k = 0; k < DIM; k++)
		{
			ah_mask_array[k] = 1;
		}
	}

	// No horizon has been found.
	ah_found = 0;

	// Set start at maximum value.
	start = ahfind_maxr;

	// Open file.
	f_H = fopen("ah.asc", "w");
}

void ah_destroy(void)
{
	// Deallocate memory.
	if (array_th)
		free(array_th);
	if (array_H)
		free(array_H);
	if (array_dH)
		free(array_dH);
	if (array_a)
		free(array_a);
	if (array_b)
		free(array_b);
	if (array_h)
		free(array_h);
	if (array_c)
		free(array_c);
	if (array_phi)
		free(array_phi);
	if (ah_mask_array)
		free(ah_mask_array);

	// Close file.
	fclose(f_H);
}

void write_ah(const double t)
{
	int k = 0;

	fprintf(f_H, "# t = %9.18E\n", t);
	for (k = 0; k < NthTotal; k++)
	{
		fprintf(f_H, "%9.18E\t%9.18E\t%9.18E\t%9.18E\t%9.18E\t%9.18E\t%9.18E\t%9.18E\t%9.18E\n", t, array_th[k], array_H[k], array_dH[k], array_a[k], array_b[k], array_h[k], array_c[k], array_phi[k]);
	}
	fprintf(f_H, "\n");
}
