#undef VERBOSE
#define MESSAGE_LEVEL 0

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h" 

#ifdef PARDISO_MAIN_FILE
// Matrix type.
MKL_INT mtype;
// Number of RHS.
MKL_INT nrhs;
// Matrix dimension.
MKL_INT n;
// Internal solver memory pointer pt.
void *pt[64];
// PARDISO control parameters.
MKL_INT iparm[64];
// Maximum number of matrices to keep in memory.
MKL_INT maxfct;
// Selected matrix to solve.
MKL_INT mnum;
// PARDISO phase.
MKL_INT phase;
// PARDISO error flag.
MKL_INT error;
// PARDISO message level;
MKL_INT msglvl;
// Dumb auxiliary variables.
double ddum;
MKL_INT idum;
// Permutation vector.
MKL_INT *perm;
// Matrix-vector multiplication type.
char *uplo;
// Low Rank parameters.
MKL_INT ndiff;
MKL_INT *diff;
#else
extern MKL_INT mtype;
extern MKL_INT nrhs;
extern MKL_INT n;
extern void *pt[64];
extern MKL_INT iparm[64];
extern MKL_INT maxfct;
extern MKL_INT mnum;
extern MKL_INT phase;
extern MKL_INT error;
extern MKL_INT msglvl;
extern double ddum;
extern MKL_INT idum;
extern MKL_INT *perm;
extern char *uplo;
extern MKL_INT ndiff;
extern MKL_INT *diff;
#endif
