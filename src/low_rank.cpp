#include "tools.h"
#include "param.h"
#include "pardiso_param.h"

// One-based indexing.
#define BASE 1

// Low rank vector print for debug.
#undef DEBUG

// Number of different elements between matrices.
int ndiff_general_elliptic(const int NrInterior, const int NzInterior, const int order)
{
	int ndiff;

	// Second order Laplacian.
	if (order == 2)
	{
		ndiff = 9 * NrInterior * NzInterior;
	}
	// Fourth order Laplacian.
	else if (order == 4)
	{
		ndiff = 17 * (NrInterior - 2) * (NzInterior - 2)
		 + (16 + 26) * (NrInterior + NzInterior - 4)
		 + (14 + 21 + 21 + 27);
	}

	return ndiff;
}

// Print low-rank vector.
void low_rank_print(const int *diff, const int ndiff, const char *fname)
{
	int k;

	FILE *fp = fopen(fname, "w");

	for (k = 0; k < 2 * ndiff + 1; k++)
		fprintf(fp, "%d\n", diff[k]);

	fclose(fp);
}

void low_rank_init(const int NrInterior, const int NzInterior, const int order)
{
	// Set number of differing elements.
	ndiff = ndiff_general_elliptic(NrInterior, NzInterior, order);

	// Allocate memory for diff array.
	diff = (int *)malloc((2 * ndiff + 1) * sizeof(int));

#ifdef VERBOSE
	printf("PARDISO LOW RANK UPDATE: Setup diff and ndiff parameters.\n");
#endif
}

void low_rank_die(void)
{
	free(diff);
}


void low_rank_set(const int NrInterior, const int NzInterior, const int order)
{
	// Auxiliary integers.
	int i, j;

	// Number of elements we have filled in.
	int offset = 0;

	// Temporary offset.
	int t_offset = 0;

	// First array of diff array is ndiff itself.
	diff[offset] = ndiff;

	// Increase offset.
	offset += 1;

	// SECOND-ORDER GENERAL ELLIPTIC EQUATION.
	if (order == 2)
	{
		// Differing elements are interior points.
		#pragma omp parallel shared(diff) private(offset, j)
		{
			#pragma omp for schedule(guided)
			for (i = 1; i < NrInterior + 1; i++)
			{
				// Each iteration of i loop will fill 2 * (9 * NzInterior) elements in diff array.
				offset = t_offset + (i - 1) * 18 * NzInterior;

				// Loop over interior points.
				for (j = 1; j < NzInterior + 1; j++)
				{
					// Row indices.
					diff[offset] 
					= diff[offset + 2] 
					= diff[offset + 4] 
					= diff[offset + 6] 
					= diff[offset + 8] 
					= diff[offset + 10] 
					= diff[offset + 12] 
					= diff[offset + 14] 
					= diff[offset + 16] = BASE + IDX(i, j);
					// Column indices.
					diff[offset + 1] = BASE + IDX(i - 1, j - 1);
					diff[offset + 3] = BASE + IDX(i - 1, j);
					diff[offset + 5] = BASE + IDX(i - 1, j + 1);
					diff[offset + 7] = BASE + IDX(i, j - 1);
					diff[offset + 9] = BASE + IDX(i, j);
					diff[offset + 11] = BASE + IDX(i, j + 1);
					diff[offset + 13] = BASE + IDX(i + 1, j - 1);
					diff[offset + 15] = BASE + IDX(i + 1, j);
					diff[offset + 17] = BASE + IDX(i + 1, j + 1);
					// Increase offset by 18.
					offset += 18;
				}
			}
		}
	}
	// FOURTH-ORDER GENERAL ELLIPTIC EQUATION.
	else if (order == 4)
	{
		// Lower-left interior corner: 2 * 14 elements.
		i = 1;
		j = 1;
		// Row indices.
		diff[offset]
		= diff[offset + 2]
		= diff[offset + 4]
		= diff[offset + 6]
		= diff[offset + 8]
		= diff[offset + 10]
		= diff[offset + 12]
		= diff[offset + 14]
		= diff[offset + 16]
		= diff[offset + 18]
		= diff[offset + 20]
		= diff[offset + 22]
		= diff[offset + 24]
		= diff[offset + 26] = BASE + IDX(i, j);
		// Column indices.
		diff[offset + 1] = BASE + IDX(0, 0);
		diff[offset + 3] = BASE + IDX(0, 1);
		diff[offset + 5] = BASE + IDX(0, 2);
		diff[offset + 7] = BASE + IDX(1, 0);
		diff[offset + 9] = BASE + IDX(1, 1);
		diff[offset + 11] = BASE + IDX(1, 2);
		diff[offset + 13] = BASE + IDX(1, 3);
		diff[offset + 15] = BASE + IDX(2, 0);
		diff[offset + 17] = BASE + IDX(2, 1);
		diff[offset + 19] = BASE + IDX(2, 2);
		diff[offset + 21] = BASE + IDX(2, 3);
		diff[offset + 23] = BASE + IDX(3, 1);
		diff[offset + 25] = BASE + IDX(3, 2);
		diff[offset + 27] = BASE + IDX(3, 3);
		// Increase offset by 28.
		offset += 28;

		// Set temporary offset.
		t_offset = offset;

		// Left interior strip: 2 * (NzInterior - 2) * 16 elements.
		#pragma omp parallel shared(diff) private(offset)
		{
			#pragma omp for schedule(guided)
			for (j = 2; j < NzInterior; j++)
			{
				// Each iteration fills 32 elements.
				offset = t_offset + (j - 2) * 32;
				// Row indices.
				diff[offset]
				= diff[offset + 2]
				= diff[offset + 4]
				= diff[offset + 6]
				= diff[offset + 8]
				= diff[offset + 10]
				= diff[offset + 12]
				= diff[offset + 14]
				= diff[offset + 16]
				= diff[offset + 18]
				= diff[offset + 20]
				= diff[offset + 22]
				= diff[offset + 24]
				= diff[offset + 26]
				= diff[offset + 28]
				= diff[offset + 30] = BASE + IDX(i, j);
				// Column indices.
				diff[offset + 1] = BASE + IDX(i - 1, j - 1);
				diff[offset + 3] = BASE + IDX(i - 1, j);
				diff[offset + 5] = BASE + IDX(i - 1, j + 1);
				diff[offset + 7] = BASE + IDX(i, j - 2);
				diff[offset + 9] = BASE + IDX(i, j - 1);
				diff[offset + 11] = BASE + IDX(i, j);
				diff[offset + 13] = BASE + IDX(i, j + 1);
				diff[offset + 15] = BASE + IDX(i, j + 2);
				diff[offset + 17] = BASE + IDX(i + 1, j - 2);
				diff[offset + 19] = BASE + IDX(i + 1, j - 1);
				diff[offset + 21] = BASE + IDX(i + 1, j);
				diff[offset + 23] = BASE + IDX(i + 1, j + 1);
				diff[offset + 25] = BASE + IDX(i + 1, j + 2);
				diff[offset + 27] = BASE + IDX(i + 2, j - 2);
				diff[offset + 29] = BASE + IDX(i + 2, j);
				diff[offset + 31] = BASE + IDX(i + 2, j + 2);
			}
		}

		// At this point we have filled:
		offset = 1 + 2 * (14 + 16 * (NzInterior - 2));

		// Top-left corner: 2 * 21 elements.
		j = NzInterior;
		// Row indices.
		diff[offset]
		= diff[offset + 2]
		= diff[offset + 4]
		= diff[offset + 6]
		= diff[offset + 8]
		= diff[offset + 10]
		= diff[offset + 12]
		= diff[offset + 14]
		= diff[offset + 16]
		= diff[offset + 18]
		= diff[offset + 20]
		= diff[offset + 22]
		= diff[offset + 24]
		= diff[offset + 26]
		= diff[offset + 28]
		= diff[offset + 30]
		= diff[offset + 32]
		= diff[offset + 34]
		= diff[offset + 36]
		= diff[offset + 38]
		= diff[offset + 40] = BASE + IDX(i, j);
		// Column indices.
		diff[offset + 1] = BASE + IDX(i - 1, j - 3);
		diff[offset + 3] = BASE + IDX(i - 1, j - 2);
		diff[offset + 5] = BASE + IDX(i - 1, j - 1);
		diff[offset + 7] = BASE + IDX(i - 1, j);
		diff[offset + 9] = BASE + IDX(i - 1, j + 1);
		diff[offset + 11] = BASE + IDX(i, j - 4);
		diff[offset + 13] = BASE + IDX(i, j - 3);
		diff[offset + 15] = BASE + IDX(i, j - 2);
		diff[offset + 17] = BASE + IDX(i, j - 1);
		diff[offset + 19] = BASE + IDX(i, j);
		diff[offset + 21] = BASE + IDX(i, j + 1);
		diff[offset + 23] = BASE + IDX(i + 1, j - 3);
		diff[offset + 25] = BASE + IDX(i + 1, j - 2);
		diff[offset + 27] = BASE + IDX(i + 1, j - 1);
		diff[offset + 29] = BASE + IDX(i + 1, j);
		diff[offset + 31] = BASE + IDX(i + 1, j + 1);
		diff[offset + 33] = BASE + IDX(i + 2, j - 3);
		diff[offset + 35] = BASE + IDX(i + 2, j - 2);
		diff[offset + 37] = BASE + IDX(i + 2, j - 1);
		diff[offset + 39] = BASE + IDX(i + 2, j);
		diff[offset + 41] = BASE + IDX(i + 2, j + 1);
		// Increase offset by 42.
		offset += 42;

		// Set temporary offset.
		t_offset = offset;

		// Interior-interior points.
		#pragma omp parallel shared(diff) private(offset, j)
		{
			#pragma omp for schedule(guided)
			for (i = 2; i < NrInterior; i++)
			{
				// Each iteration of i loop will fill 2 * (16 + 17 * (NzInterior - 2) + 26) elements in diff array.
				offset = t_offset + (i - 2) * (16 + 34 * NzInterior);

				// Lower strip: 2 * 16 elements.
				j = 1;
				// Row indices.
				diff[offset]
				= diff[offset + 2]
				= diff[offset + 4]
				= diff[offset + 6]
				= diff[offset + 8]
				= diff[offset + 10]
				= diff[offset + 12]
				= diff[offset + 14]
				= diff[offset + 16]
				= diff[offset + 18]
				= diff[offset + 20]
				= diff[offset + 22]
				= diff[offset + 24]
				= diff[offset + 26]
				= diff[offset + 28]
				= diff[offset + 30] = BASE + IDX(i, j);
				// Column indices.
				diff[offset + 1] = BASE + IDX(i - 2, j);
				diff[offset + 3] = BASE + IDX(i - 2, j + 1);
				diff[offset + 5] = BASE + IDX(i - 2, j + 2);
				diff[offset + 7] = BASE + IDX(i - 1, j - 1);
				diff[offset + 9] = BASE + IDX(i - 1, j);
				diff[offset + 11] = BASE + IDX(i - 1, j + 1);
				diff[offset + 13] = BASE + IDX(i, j - 1);
				diff[offset + 15] = BASE + IDX(i, j);
				diff[offset + 17] = BASE + IDX(i, j + 1);
				diff[offset + 19] = BASE + IDX(i, j + 2);
				diff[offset + 21] = BASE + IDX(i + 1, j - 1);
				diff[offset + 23] = BASE + IDX(i + 1, j);
				diff[offset + 25] = BASE + IDX(i + 1, j + 1);
				diff[offset + 27] = BASE + IDX(i + 2, j);
				diff[offset + 29] = BASE + IDX(i + 2, j + 1);
				diff[offset + 31] = BASE + IDX(i + 2, j + 2);
				// Increase offset by 32.
				offset += 32;

				// Loop over interior-interior points: 2 * (NzInterior - 2) * 17 elements.
				for (j = 2; j < NzInterior; j++)
				{
					// Row indices.
					diff[offset]
					= diff[offset + 2]
					= diff[offset + 4]
					= diff[offset + 6]
					= diff[offset + 8]
					= diff[offset + 10]
					= diff[offset + 12]
					= diff[offset + 14]
					= diff[offset + 16]
					= diff[offset + 18]
					= diff[offset + 20]
					= diff[offset + 22]
					= diff[offset + 24]
					= diff[offset + 26]
					= diff[offset + 28]
					= diff[offset + 30]
					= diff[offset + 32] = BASE + IDX(i, j);
					// Column indices.
					diff[offset + 1] = BASE + IDX(i - 2, j - 2);
					diff[offset + 3] = BASE + IDX(i - 2, j);
					diff[offset + 5] = BASE + IDX(i - 2, j + 2);
					diff[offset + 7] = BASE + IDX(i - 1, j - 1);
					diff[offset + 9] = BASE + IDX(i - 1, j);
					diff[offset + 11] = BASE + IDX(i - 1, j + 1);
					diff[offset + 13] = BASE + IDX(i, j - 2);
					diff[offset + 15] = BASE + IDX(i, j - 1);
					diff[offset + 17] = BASE + IDX(i, j);
					diff[offset + 19] = BASE + IDX(i, j + 1);
					diff[offset + 21] = BASE + IDX(i, j + 2);
					diff[offset + 23] = BASE + IDX(i + 1, j - 1);
					diff[offset + 25] = BASE + IDX(i + 1, j);
					diff[offset + 27] = BASE + IDX(i + 1, j + 1);
					diff[offset + 29] = BASE + IDX(i + 2, j - 2);
					diff[offset + 31] = BASE + IDX(i + 2, j);
					diff[offset + 33] = BASE + IDX(i + 2, j + 2);
					// Increase offset by 34.
					offset += 34;
				}

				// Top strip: 2 * 26 elements.
				j = NzInterior;
				// Row indices.
				diff[offset]
				= diff[offset + 2]
				= diff[offset + 4]
				= diff[offset + 6]
				= diff[offset + 8]
				= diff[offset + 10]
				= diff[offset + 12]
				= diff[offset + 14]
				= diff[offset + 16]
				= diff[offset + 18]
				= diff[offset + 20]
				= diff[offset + 22]
				= diff[offset + 24]
				= diff[offset + 26]
				= diff[offset + 28]
				= diff[offset + 30]
				= diff[offset + 32]
				= diff[offset + 34]
				= diff[offset + 36]
				= diff[offset + 38]
				= diff[offset + 40]
				= diff[offset + 42]
				= diff[offset + 44]
				= diff[offset + 46]
				= diff[offset + 48]
				= diff[offset + 50] = BASE + IDX(i, j);
				// Column indices.
				diff[offset + 1] = BASE + IDX(i - 2, j - 3);
				diff[offset + 3] = BASE + IDX(i - 2, j - 2);
				diff[offset + 5] = BASE + IDX(i - 2, j - 1);
				diff[offset + 7] = BASE + IDX(i - 2, j);
				diff[offset + 9] = BASE + IDX(i - 2, j + 1);
				diff[offset + 11] = BASE + IDX(i - 1, j - 3);
				diff[offset + 13] = BASE + IDX(i - 1, j - 2);
				diff[offset + 15] = BASE + IDX(i - 1, j - 1);
				diff[offset + 17] = BASE + IDX(i - 1, j);
				diff[offset + 19] = BASE + IDX(i - 1, j + 1);
				diff[offset + 21] = BASE + IDX(i, j - 4);
				diff[offset + 23] = BASE + IDX(i, j - 3);
				diff[offset + 25] = BASE + IDX(i, j - 2);
				diff[offset + 27] = BASE + IDX(i, j - 1);
				diff[offset + 29] = BASE + IDX(i, j);
				diff[offset + 31] = BASE + IDX(i, j + 1);
				diff[offset + 33] = BASE + IDX(i + 1, j - 3);
				diff[offset + 35] = BASE + IDX(i + 1, j - 2);
				diff[offset + 37] = BASE + IDX(i + 1, j - 1);
				diff[offset + 39] = BASE + IDX(i + 1, j);
				diff[offset + 41] = BASE + IDX(i + 1, j + 1);
				diff[offset + 43] = BASE + IDX(i + 2, j - 3);
				diff[offset + 45] = BASE + IDX(i + 2, j - 2);
				diff[offset + 47] = BASE + IDX(i + 2, j - 1);
				diff[offset + 49] = BASE + IDX(i + 2, j);
				diff[offset + 51] = BASE + IDX(i + 2, j + 1);
				// offset += 52;
			}
		}

		// At this point we have filled:
		offset = 1 + 2 * (14 + 16 * (NzInterior - 2) + 21 + 16 * (NrInterior - 2) + 17 * (NrInterior - 2) * (NzInterior - 2) + 26 * (NrInterior - 2));

		// Interior bottom-right corner: 2 * 21 elements.	
		i = NrInterior;
		j = 1;
		// Row indices.
		diff[offset]
		= diff[offset + 2]
		= diff[offset + 4]
		= diff[offset + 6]
		= diff[offset + 8]
		= diff[offset + 10]
		= diff[offset + 12]
		= diff[offset + 14]
		= diff[offset + 16]
		= diff[offset + 18]
		= diff[offset + 20]
		= diff[offset + 22]
		= diff[offset + 24]
		= diff[offset + 26]
		= diff[offset + 28]
		= diff[offset + 30]
		= diff[offset + 32]
		= diff[offset + 34]
		= diff[offset + 36]
		= diff[offset + 38]
		= diff[offset + 40] = BASE + IDX(i, j);
		// Column indices.
		diff[offset + 1] = BASE + IDX(i - 4, j);
		diff[offset + 3] = BASE + IDX(i - 3, j - 1);
		diff[offset + 5] = BASE + IDX(i - 3, j);
		diff[offset + 7] = BASE + IDX(i - 3, j + 1);
		diff[offset + 9] = BASE + IDX(i - 3, j + 2);
		diff[offset + 11] = BASE + IDX(i - 2, j - 1);
		diff[offset + 13] = BASE + IDX(i - 2, j);
		diff[offset + 15] = BASE + IDX(i - 2, j + 1);
		diff[offset + 17] = BASE + IDX(i - 2, j + 2);
		diff[offset + 19] = BASE + IDX(i - 1, j - 1);
		diff[offset + 21] = BASE + IDX(i - 1, j);
		diff[offset + 23] = BASE + IDX(i - 1, j + 1);
		diff[offset + 25] = BASE + IDX(i - 1, j + 2);
		diff[offset + 27] = BASE + IDX(i, j - 1);
		diff[offset + 29] = BASE + IDX(i, j);
		diff[offset + 31] = BASE + IDX(i, j + 1);
		diff[offset + 33] = BASE + IDX(i, j + 2);
		diff[offset + 35] = BASE + IDX(i + 1, j - 1);
		diff[offset + 37] = BASE + IDX(i + 1, j);
		diff[offset + 39] = BASE + IDX(i + 1, j + 1);
		diff[offset + 41] = BASE + IDX(i + 1, j + 2);
		// Increase offset by 42.
		offset += 42;

		// Set temporary offset.
		t_offset = offset;

		// Right strip: 2 * (NzInterior - 2) * 26 elements.
		#pragma omp parallel shared(diff) private(offset)
		{
			#pragma omp for schedule(guided)
			for (j = 2; j < NzInterior; j++)
			{
				// Each iteration fills 52 elements.
				offset = t_offset + (j - 2) * 52;
				// Row indices.
				diff[offset]
				= diff[offset + 2]
				= diff[offset + 4]
				= diff[offset + 6]
				= diff[offset + 8]
				= diff[offset + 10]
				= diff[offset + 12]
				= diff[offset + 14]
				= diff[offset + 16]
				= diff[offset + 18]
				= diff[offset + 20]
				= diff[offset + 22]
				= diff[offset + 24]
				= diff[offset + 26]
				= diff[offset + 28]
				= diff[offset + 30]
				= diff[offset + 32]
				= diff[offset + 34]
				= diff[offset + 36]
				= diff[offset + 38]
				= diff[offset + 40]
				= diff[offset + 42]
				= diff[offset + 44]
				= diff[offset + 46]
				= diff[offset + 48]
				= diff[offset + 50] = BASE + IDX(i, j);
				// Column indices.
				diff[offset + 1] = BASE + IDX(i - 4, j);
				diff[offset + 3] = BASE + IDX(i - 3, j - 2);
				diff[offset + 5] = BASE + IDX(i - 3, j - 1);
				diff[offset + 7] = BASE + IDX(i - 3, j);
				diff[offset + 9] = BASE + IDX(i - 3, j + 1);
				diff[offset + 11] = BASE + IDX(i - 3, j + 2);
				diff[offset + 13] = BASE + IDX(i - 2, j - 2);
				diff[offset + 15] = BASE + IDX(i - 2, j - 1);
				diff[offset + 17] = BASE + IDX(i - 2, j);
				diff[offset + 19] = BASE + IDX(i - 2, j + 1);
				diff[offset + 21] = BASE + IDX(i - 2, j + 2);
				diff[offset + 23] = BASE + IDX(i - 1, j - 2);
				diff[offset + 25] = BASE + IDX(i - 1, j - 1);
				diff[offset + 27] = BASE + IDX(i - 1, j);
				diff[offset + 29] = BASE + IDX(i - 1, j + 1);
				diff[offset + 31] = BASE + IDX(i - 1, j + 2);
				diff[offset + 33] = BASE + IDX(i, j - 2);
				diff[offset + 35] = BASE + IDX(i, j - 1);
				diff[offset + 37] = BASE + IDX(i, j);
				diff[offset + 39] = BASE + IDX(i, j + 1);
				diff[offset + 41] = BASE + IDX(i, j + 2);
				diff[offset + 43] = BASE + IDX(i + 1, j - 2);
				diff[offset + 45] = BASE + IDX(i + 1, j - 1);
				diff[offset + 47] = BASE + IDX(i + 1, j);
				diff[offset + 49] = BASE + IDX(i + 1, j + 1);
				diff[offset + 51] = BASE + IDX(i + 1, j + 2);
			}
		}

		// At this point we have filled:
		offset = 1 + 2 * (14 + 16 * (NzInterior - 2) + 21 + 16 * (NrInterior - 2) + 17 * (NrInterior - 2) * (NzInterior - 2) + 26 * (NrInterior - 2) + 21 + 26 * (NzInterior - 2));

		// Interior top-right corner: 2 * 27 elements.
		j = NzInterior;
		// Row indices.
		diff[offset]
		= diff[offset + 2]
		= diff[offset + 4]
		= diff[offset + 6]
		= diff[offset + 8]
		= diff[offset + 10]
		= diff[offset + 12]
		= diff[offset + 14]
		= diff[offset + 16]
		= diff[offset + 18]
		= diff[offset + 20]
		= diff[offset + 22]
		= diff[offset + 24]
		= diff[offset + 26]
		= diff[offset + 28]
		= diff[offset + 30]
		= diff[offset + 32]
		= diff[offset + 34]
		= diff[offset + 36]
		= diff[offset + 38]
		= diff[offset + 40]
		= diff[offset + 42]
		= diff[offset + 44]
		= diff[offset + 46]
		= diff[offset + 48]
		= diff[offset + 50]
		= diff[offset + 52] = BASE + IDX(i, j);
		// Columns.
		diff[offset + 1] = BASE + IDX(i - 4, j);
		diff[offset + 3] = BASE + IDX(i - 3, j - 3);
		diff[offset + 5] = BASE + IDX(i - 3, j - 2);
		diff[offset + 7] = BASE + IDX(i - 3, j - 1);
		diff[offset + 9] = BASE + IDX(i - 3, j);
		diff[offset + 11] = BASE + IDX(i - 3, j + 1);
		diff[offset + 13] = BASE + IDX(i - 2, j - 3);
		diff[offset + 15] = BASE + IDX(i - 2, j - 2);
		diff[offset + 17] = BASE + IDX(i - 2, j - 1);
		diff[offset + 19] = BASE + IDX(i - 2, j);
		diff[offset + 21] = BASE + IDX(i - 2, j + 1);
		diff[offset + 23] = BASE + IDX(i - 1, j - 3);
		diff[offset + 25] = BASE + IDX(i - 1, j - 2);
		diff[offset + 27] = BASE + IDX(i - 1, j - 1);
		diff[offset + 29] = BASE + IDX(i - 1, j);
		diff[offset + 31] = BASE + IDX(i - 1, j + 1);
		diff[offset + 33] = BASE + IDX(i, j - 4);
		diff[offset + 35] = BASE + IDX(i, j - 3);
		diff[offset + 37] = BASE + IDX(i, j - 2);
		diff[offset + 39] = BASE + IDX(i, j - 1);
		diff[offset + 41] = BASE + IDX(i, j);
		diff[offset + 43] = BASE + IDX(i, j + 1);
		diff[offset + 45] = BASE + IDX(i + 1, j - 3);
		diff[offset + 47] = BASE + IDX(i + 1, j - 2);
		diff[offset + 49] = BASE + IDX(i + 1, j - 1);
		diff[offset + 51] = BASE + IDX(i + 1, j);
		diff[offset + 53] = BASE + IDX(i + 1, j + 1);
		// offset += 54;
	}

#ifdef DEBUG
	low_rank_print(diff, ndiff, "diff.asc");
#endif

	// All done.
	return;
}
