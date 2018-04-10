#include "tools.h"

#include "param.h"
#include "arrays.h"

#include "radbound.h"

void boundary(void)
{
	/* RADIATIVE BOUNDARIES */
	double speed = 1.;
	double var0 = 0.;
	double var1 = 1.;

	if (strcmp(boundtype, "radiative") == 0)
	{
		// Radiative boundary for traceless extrinsic curvature.
		radbound(speed, var0, A_a, Dr_A_a, Dz_A_a, sA_a);
		radbound(speed, var0, A_b, Dr_A_b, Dz_A_b, sA_b);
		radbound(speed, var0, A_h, Dr_A_h, Dz_A_h, sA_h);
		radbound(speed, var0, A_c, Dr_A_c, Dz_A_c, sA_c);
		radbound(speed, var0, A_lambda, Dr_A_lambda, Dz_A_lambda, sA_lambda);

#ifdef SHIFT
		// Shift terms.
		if (!(strcmp(shift, "none") == 0))
		{
			radbound(speed, var0, Deltar, Dr_Deltar, Dz_Deltar, sDeltar);
			radbound(speed, var0, Deltaz, Dr_Deltaz, Dz_Deltaz, sDeltaz);

			radbound(1.0, 0.0, dtbeta_r, Dr_dtbeta_r, Dz_dtbeta_r, sdtbeta_r);
			radbound(1.0, 0.0, dtbeta_z, Dr_dtbeta_z, Dz_dtbeta_z, sdtbeta_z);
		}
#endif

		// Radiative boundary for trK.
		if (!(strcmp(slicing, "maximal") == 0))
		{
			speed = sqrt(gauge_f);
			radbound(speed, var0, K, Dr_K, Dz_K, sK);
		}
	}

#ifdef DISS
	/* DISSIPATION */
	if (diss != 0.0)
	{
		dissipation(a, sa, 1, 1);
		dissipation(b, sb, 1, 1);
		dissipation(h, sh, 1, 1);
		dissipation(c, sc, 1, -1);
		dissipation(lambda, slambda, 1, 1);
		if (!(strcmp(slicing, "maximal") == 0))
		{
			dissipation(phi, sphi, 1, 1);
			dissipation(K, sK, 1, 1);
		}
		dissipation(A_a, sA_a, 1, 1);
		dissipation(A_b, sA_b, 1, 1);
		dissipation(A_h, sA_h, 1, 1);
		dissipation(A_c, sA_c, 1, -1);
		dissipation(A_lambda, sA_lambda, 1, 1);
		dissipation(Deltar, sDeltar, -1, 1);
		dissipation(Deltaz, sDeltaz, 1, -1);
	}
#endif

	return;
}