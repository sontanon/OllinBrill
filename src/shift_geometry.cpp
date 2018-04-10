#include "tools.h"

#include "param.h"
#include "arrays.h"

#include "derivatives.h"
#include "shift_derivadvect.h"
#include "shift_divergence.h"

void shift_geometry(void)
{
	// Auxiliary integers.
	int i, j, k;

	// The determinant is one, so its derivatives are zero.
	// This quantity only evolves for a non-zero eulerian shift.
	if ((strcmp(bssn_flavor, "eulerian") == 0))
	{
		diff1r(Dr_hdet, hdet, 1);
		diff1z(Dz_hdet, hdet, 1);
	}

	// Shift derivatives.
	diff1r(Dr_beta_r, beta_r, -1);
	diff1z(Dz_beta_r, beta_r, 1);
	diff2r(Drr_beta_r, beta_r, -1);
	diff2z(Dzz_beta_r, beta_r, 1);
	diff2rz(Drz_beta_r, beta_r, -1, 1);
	diff1r(Dr_beta_z, beta_z, 1);
	diff1z(Dz_beta_z, beta_z, -1);
	diff2r(Drr_beta_z, beta_z, 1);
	diff2z(Dzz_beta_z, beta_z, -1);
	diff2rz(Drz_beta_z, beta_z, 1, -1);
	// Regularization.
	#pragma omp parallel shared(auxarray)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			auxarray[k] = beta_r[k] / r[k];
		}
	}
	diff1r(DD_beta_rr, auxarray, 1);

	// Time derivatives.
	diff1r(Dr_dtbeta_r, dtbeta_r, -1);
	diff1z(Dz_dtbeta_r, dtbeta_r, 1);
	diff1r(Dr_dtbeta_z, dtbeta_z, 1);
	diff1z(Dz_dtbeta_z, dtbeta_z, -1);

	// Divergence.
	#pragma omp parallel shared(div_beta)
	{
		#pragma omp for schedule(guided)
		for (k = 0; k < DIM; k++)
		{
			//beta_div(div_beta + k, r[k], beta_r[k], beta_z[k], Dr_beta_r[k], Dz_beta_z[k], hdet[k], Dr_hdet[k], Dz_hdet[k]);
			beta_div(div_beta + k, r[k], beta_r[k], beta_z[k], Dr_beta_r[k], Dz_beta_z[k], 1.0, 0.0, 0.0);
		}
	}
	diff1r(Dr_div_beta, div_beta, 1);
	diff1z(Dz_div_beta, div_beta, 1);

	// Advective derivatives.
	shift_diffadvr();
	shift_diffadvz();

	return;
}