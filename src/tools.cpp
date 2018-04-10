// Windows architecture.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <omp.h>
#include <string.h>

// LIBCONFIG.
#include <libconfig.h>

// Intel MKL headers.
#include "mkl.h"

// Global variables.
#include "param.h"

// CSR Matrix index type.
#define BASE 1
//#define ALIGN 128

// Types.
typedef struct csr_matrices
{
	double *a;
	double *b;
	int *ia;
	int *ja;
	int nrows;
	int ncols;
	int nnz;

} csr_matrix;

// Forward declarations.
int IDX(const int i, const int j);

// Indexing macro.
int IDX(const int i, const int j)
{
	return i * NzTotal + j;
}

void writeSingleFile(const double *u, const char *fname);

// CSR matrix print.
void csr_print(csr_matrix A, const char *vA, const char *iA, const char *jA)
{
	FILE *fvA = fopen(vA, "w");
	FILE *fiA = fopen(iA, "w");
	FILE *fjA = fopen(jA, "w");

	int k = 0;

	for (k = 0; k < A.nnz; k++)
	{
		fprintf(fvA, "%9.18E\n", A.a[k]);
		fprintf(fjA, "%d\n", A.ja[k]);
	}

	for (k = 0; k < A.nrows + 1; k++)
		fprintf(fiA, "%d\n", A.ia[k]);

	fclose(fvA);
	fclose(fiA);
	fclose(fjA);
}

void writeSingleFile(const double *u, const char *fname)
{
	int i, j;

	FILE *fp = fopen(fname, "w");

	for (i = 0; i < NrTotal; i++)
	{
		for (j = 0; j < NzTotal; j++)
		{
			fprintf(fp, (j < NzTotal - 1) ? "%9.18E\t" : "%9.18E\n", u[IDX(i, j)]);
		}
	}

	fclose(fp);
}

// Construct CSR matrix.
void csr_allocate(csr_matrix *A, const int nrows, const int ncols,
	const int nnz)
{
	A->nrows = nrows;
	A->ncols = ncols;
	A->nnz = nnz;
	A->a = (double *)malloc(sizeof(double) * nnz);
	A->ja = (int *)malloc(sizeof(int) * nnz);
	A->ia = (int *)malloc(sizeof(int) * (nrows + 1));
}

void csr_deallocate(csr_matrix *A)
{
	A->nrows = 0;
	A->ncols = 0;
	A->nnz = 0;
	free(A->a);
	free(A->ia);
	free(A->ja);
}
