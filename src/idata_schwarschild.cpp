#include "tools.h"

#include "param.h"
#include "arrays.h"

#include "derivatives.h"

void idata_schwarschild(void)
{
	// Auxiliary variables.
	int i, j, k;
	double aux_r, aux_z, aux_rr, aux_1;

#pragma omp parallel private(j, aux_r, aux_z, aux_rr, aux_1)\
	shared(phi, psi, psi4)
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

				// Radial coordinate.
				aux_rr = rr[IDX(i, j)];

				// Conformal factor.
				aux_1 = (1.0 + 0.5 * schwar_mass / aux_rr);

				// Write to arrays.
				psi[IDX(i, j)] = aux_1;
				psi4[IDX(i, j)] = aux_1 * aux_1 * aux_1 * aux_1;
				phi[IDX(i, j)] = log(ABS(aux_1));
			}
		}
	}

	// Calculate conformal factor derivatives.
	diff1r(Dr_phi, phi, 1);
	diff1z(Dz_phi, phi, 1);
	// Regularization.
	vdDiv(DIM, Dr_phi, r, auxarray);
	diff1r(DDphi_r, auxarray, 1);
	// Second derivatives.
	diff2r(Drr_phi, phi, 1);
	diff2z(Dzz_phi, phi, 1);
	diff2rz(Drz_phi, phi, 1, 1);

	// Print out Brill parameters.
	printf("***                                ***\n");
	printf("***   Generated initial data for   ***\n");
	printf("***    Schwarschild black hole:    ***\n");
	printf("***                                ***\n");
	printf("***   schwar_mass =  %-10.3E     ***\n", schwar_mass);
	printf("***                                ***\n");
	printf("**************************************\n");
}
