#include "tools.h"

#include "param.h"
#include "arrays.h"
#include "derivatives.h"
#include "flat_laplacian.h"
#include "pardiso_init.h"
#include "pardiso_die.h"

// Define to print out elliptic equation.
#undef DEBUG
#define DEBUG0

void idata_brillwave(void)
{
	// Auxiliary variables.
	int i, j, k;
	double aux_r, aux_z, aux1, aux2;

	#pragma omp parallel private(j, aux_r, aux_z, aux1, aux2) \
	shared(auxarray, auxarray2, auxarray3, auxarray4, a, b, h,\
	lambda, Deltar, Deltaz, A_a, A_b, A_h, A_c, A_lambda) 	
	{
		#pragma omp for schedule(guided)
		for (i = 0; i < NrTotal; i++)
		{
			// R coordinate.
			aux_r = r[IDX(i, 0)];

			// Loop over Z.
			for (j = 0; j < NzTotal; j++)
			{
				// Z coordinate.
				aux_z = z[IDX(i, j)];

				// Calculate q function omiiting r**2 factor.
				aux1 = brill_a0 * exp(-pow(aux_r / brill_sr0, 2))
					* exp(-pow(aux_z / brill_sz0, 2));

				// Write q function to auxarray.
				auxarray[IDX(i, j)] = aux1 * pow(aux_r, 2);

				// Linear source for Laplacian, write to auxarray2.
				auxarray2[IDX(i, j)] = aux1 * (0.5 + pow(aux_r, 2)
					* (-2.5 / pow(brill_sr0, 2) - 0.5 / pow(brill_sz0, 2)
						+ pow(aux_r, 2) / pow(brill_sr0, 4)
						+ pow(aux_z, 2) / pow(brill_sz0, 4)));

				// Third array serves as RHS.
				auxarray3[IDX(i, j)] = 0.0;

				// Fourth array serves as residual.
				auxarray4[IDX(i, j)] = 0.0;

				// Metric coefficients, inserting ommited factor.
				// A = B = exp(2q/3).
				aux2 = exp(2.0 * aux1 * pow(aux_r, 2) / 3.0);
				a[IDX(i, j)] = aux2;
				b[IDX(i, j)] = aux2;
				// H = exp(-4q/3).
				h[IDX(i, j)] = 1.0 / pow(aux2, 2);
				// C is set to zero by Minkowski.

				// Lambda regularization.
				lambda[IDX(i, j)] = (aux2 - 1.0 / pow(aux2, 2)) / pow(aux_r, 2);

				// Delta vector.
				Deltar[IDX(i, j)] = (1.0 / aux2) *
					(((pow(aux2, 3) - 1.0) / aux_r)
						+ (4.0 / 3.0) * aux1 * aux_r
						* (1.0 - pow(aux_r / brill_sr0, 2)));
				Deltaz[IDX(i, j)] = -(4.0 / (3.0 * aux2))
					* aux_z * aux1 * pow(aux_r / brill_sz0, 2);
			}
		}
	}

	printf("**************************************\n");
	printf("***                                ***\n");
	printf("***        Calculating Brill       ***\n");
	printf("***         initial data...        ***\n");
	printf("***                                ***\n");

	// Initialize solver.
	pardiso_init(NrInterior, NzInterior);

	// Call elliptic solver.
	flat_laplacian(psi, auxarray3, auxarray4, auxarray2, 1.0, nrobin, 1, 1, 0, 0, 0);

	// Delete solver memory.
	pardiso_die();

#ifdef DEBUG0
	writeSingleFile(auxarray, "q.asc");
	writeSingleFile(psi, "psi_bar.asc");
#endif

#ifdef DEBUG
	// Write psi.
	writeSingleFile(psi, "psi_bar.asc");

	// Write q function.
	writeSingleFile(auxarray, "q.asc");

	// Write linear source.
	writeSingleFile(auxarray2, "s.asc");

	// Print out derivatives.
	// Used so we can check the full equation:
	// Drr_psi + Dzz_psi + (1/r)Dr_psi + s * psi = 0.
	diff1r(auxarray2, psi, 1);
	diff2r(auxarray3, psi, 1);
	writeSingleFile(auxarray2, "Dr_psi_bar.asc");
	writeSingleFile(auxarray3, "Drr_psi_bar.asc");

	diff2z(auxarray3, psi, 1);
	writeSingleFile(auxarray3, "Dzz_psi_bar.asc");
#endif

	// Correct psi by exp(q/3) and calculate phi.
	#pragma omp parallel shared(psi)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			psi[k] *= exp(auxarray[k] / 3.0);
		}
	}

	// Calculate conformal factor phi.
	#pragma omp parallel shared(phi, psi4) private(aux1)
	{
		#pragma omp parallel for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			aux1 = psi[k];
			phi[k] = log(aux1);
			psi4[k] = pow(aux1, 4);
		}
	}

	// Calculate conformal factor derivatives.
	diff1r(Dr_phi, phi, 1);
	diff1z(Dz_phi, phi, 1);
	// Regularization.
	#pragma omp parallel shared(auxarray)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			auxarray[k] = Dr_phi[k] / r[k];
		}
	}
	diff1r(DDphi_r, auxarray, 1);
	// Second derivatives.
	diff2r(Drr_phi, phi, 1);
	diff2z(Dzz_phi, phi, 1);
	diff2rz(Drz_phi, phi, 1, 1);

	// Print out Brill parameters.
	printf("***                                ***\n");
	printf("***   Generated initial data for   ***\n");
	printf("***   Brill wave with parameters   ***\n");
	printf("***                                ***\n");
	printf("***     brill_a0  = %-10.3E     ***\n", brill_a0);
	printf("***     brill_sr0 = %-10.3E     ***\n", brill_sr0);
	printf("***     brill_sz0 = %-10.3E     ***\n", brill_sz0);
	printf("***                                ***\n");
	printf("**************************************\n");

	return;
}
