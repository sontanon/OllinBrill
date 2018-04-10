// AH Finder.
#define AH

// SHIFT.
#define SHIFT

// Architecture
#define WIN

// Standard headers.
#include <stdio.h>
#include <stdlib.h>
#ifdef WIN
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <time.h>
#include <assert.h>
#include <omp.h>
#include <string.h>

// LIBCONFIG
#include <libconfig.h>

// System headers.
#include <sys/types.h>
#include <sys/stat.h>
#ifdef WIN
#include <io.h>
#include <direct.h>
#else
#include <unistd.h>
#endif

// Intel MKL headers.
#include "mkl.h"

// MIN/MAX macros.
#define MIN(X, Y) ((X) < (Y)) ? (X) : (Y)
#define MAX(X, Y) ((X) > (Y)) ? (X) : (Y)

// ABS macro.
#define ABS(X) ((X) < 0) ? -(X) : (X)

//#define ALIGN 128
#define BASE 1

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
void csr_print(csr_matrix A, const char *vA, const char *iA, const char *jA);
void writeSingleFile(const double *u, const char *fname);
void csr_allocate(csr_matrix *A, const int nrows, const int ncols, const int nnz);
void csr_deallocate(csr_matrix *A);
