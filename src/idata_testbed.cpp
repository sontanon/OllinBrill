#include "tools.h"

#include "param.h"
#include "arrays.h"
#include "derivatives.h"

void idata_testbed(void)
{
	int i, j, k;
	double aux_r, aux_z;

	#pragma omp parallel private(j, aux_r, aux_z)\
	shared(a, b, h, c, lambda, Deltar, Deltaz,\
	A_a, A_b, A_h, A_c, A_lambda, phi, psi, psi4)
	{
		#pragma omp for schedule(guided)
		for (i = 0; i < NrTotal; i++)
		{
			aux_r = r[IDX(i, 0)];

			for (j = 0; j < NzTotal; j++)
			{
				aux_z = z[IDX(i, j)];

				// Metric functions.
				a[IDX(i, j)] = (exp((2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3.));
				b[IDX(i, j)] = (exp((2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3.)*
					(1 + exp(-2 * (aux_r*aux_r) -
					(4 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3. -
						2 * (aux_z*aux_z))*(aux_r*aux_r)*(aux_z*aux_z)));
				h[IDX(i, j)] = (exp((-4 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3.));
				c[IDX(i, j)] = exp(-(aux_r*aux_r) - aux_z * aux_z)*aux_z;
				lambda[IDX(i, j)] = (-exp((-4 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3.) + exp((2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3.)) / (aux_r*aux_r);

				// Delta vector.
				Deltar[IDX(i, j)] = ((exp((-3 - 2 * exp(-(aux_r*aux_r) - aux_z * aux_z))*(aux_r*aux_r) -
						3 * (aux_z*aux_z))*aux_r*(-4 *
						exp((2 +
						(4 * exp(-(aux_r*aux_r) - aux_z * aux_z)) / 3.)*(aux_r*aux_r)\
							+ 2 * (aux_z*aux_z))*(-1 + aux_r * aux_r) +
						8 * exp(aux_r*aux_r +
						(2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3. +
							aux_z * aux_z)*(aux_r*aux_r)*(aux_z*aux_z) -
						12 * (aux_r*aux_r)*(-1 + aux_r * aux_r)*(aux_z*aux_z) +
						3 * exp(aux_r*aux_r + aux_z * aux_z)*(-3 + 4 * (aux_r*aux_r))*
						(aux_z*aux_z) + exp(
						(2 + (2 * exp(-(aux_r*aux_r) - aux_z * aux_z)) / 3.)*
							(aux_r*aux_r) + 2 * (aux_z*aux_z))*(3 - 6 * (aux_z*aux_z)))) / 3.
					+ ((-1 + exp(2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r))) / aux_r) /
					(exp((2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3.)));
				Deltaz[IDX(i, j)] = ((-2 * exp((-2 - (4 * exp(-(aux_r*aux_r) - aux_z * aux_z)) / 3.)*
					(aux_r*aux_r) - 2 * (aux_z*aux_z))*
					(2 * exp(aux_r*aux_r +
					(2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3. +
						aux_z * aux_z)*(aux_r*aux_r) +
						3 * exp(aux_r*aux_r + aux_z * aux_z)*(-1 + aux_r * aux_r) -
						4 * (aux_r*aux_r)*(-1 + aux_r * aux_r))*aux_z) / 3.);

				// Curvature.
				A_a[IDX(i, j)] = (1 - (exp((2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3.)*
					(2 / exp((2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3.) + exp((4 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3.) +
					(aux_r*aux_r*(aux_z*aux_z)) / exp(2 * ((1 + exp(-(aux_r*aux_r) - aux_z * aux_z))*(aux_r*aux_r) + aux_z * aux_z)))) / 3.);
				A_b[IDX(i, j)] = ((1 - (exp((2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3.)*
					(1 + exp(-2 * (aux_r*aux_r) - (4 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3. - 2 * (aux_z*aux_z))*(aux_r*aux_r)*(aux_z*aux_z))*
					(2 / exp((2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3.) + exp((4 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3.) +
					(aux_r*aux_r*(aux_z*aux_z)) / exp(2 * ((1 + exp(-(aux_r*aux_r) - aux_z * aux_z))*(aux_r*aux_r) + aux_z * aux_z)))) / 3.));
				A_h[IDX(i, j)] = ((2 - 2 / exp(2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) -
					exp((-2 - (10 * exp(-(aux_r*aux_r) - aux_z * aux_z)) / 3.)*(aux_r*aux_r) - 2 * (aux_z*aux_z))*(aux_r*aux_r)*(aux_z*aux_z)) / 3.);
				A_c[IDX(i, j)] = (-(exp(-(aux_r*aux_r) - aux_z * aux_z)*aux_z*(2 / exp((2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3.) +
					exp((4 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r)) / 3.) +
					(aux_r*aux_r*(aux_z*aux_z)) / exp(2 * ((1 + exp(-(aux_r*aux_r) - aux_z * aux_z))*(aux_r*aux_r) + aux_z * aux_z)))) / 3.);
				A_lambda[IDX(i, j)] = (-(exp(2 * (aux_z*aux_z)) + 2 * exp(-2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r) + 2 * (aux_z*aux_z)) +
					exp((-2 - (10 * exp(-(aux_r*aux_r) - aux_z * aux_z)) / 3.)*(aux_r*aux_r))*(aux_r*aux_r)*(aux_z*aux_z))*((-1 + exp(2 * exp(-(aux_r*aux_r) - aux_z * aux_z)*(aux_r*aux_r))) / (aux_r*aux_r)) / (3.*exp(2 * (aux_z*aux_z))));

				// Conformal factor.
				phi[IDX(i, j)] = 0.0;
				psi[IDX(i, j)] = 1.0;
				psi4[IDX(i, j)] = 1.0;
			}
		}
	}

	return;
}