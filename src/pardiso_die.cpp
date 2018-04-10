#include "tools.h"
#include "pardiso_param.h"

void pardiso_die(void)
{
	// Termination and release of memory.
	phase = -1;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, &idum, &idum, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

	// Delete permutation vector.
	free(perm);

	return;
}
