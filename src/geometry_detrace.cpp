#include "tools.h"

#include "param.h"
#include "arrays.h"

void detrace_A(void)
{
	int k;
	double aux;
	const double third = 1. / 3.;

	#pragma omp parallel shared(A_a, A_b, A_h, A_c, A_lambda) private(aux)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			// Take trace.
			aux = g_a[k] * A_a[k] + g_b[k] * A_b[k] + g_h[k] * A_h[k] + 2. * r[k] * r[k] * g_c[k] * A_c[k];

			// Correct.
			A_a[k] -= third * aux * a[k];
			A_b[k] -= third * aux * b[k];
			A_h[k] -= third * aux * h[k];
			A_c[k] -= third * aux * c[k];
			A_lambda[k] -= third * aux * lambda[k];
		}
	}

	return;
}
