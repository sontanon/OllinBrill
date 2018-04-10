#include "tools.h"

#include "param.h"
#include "arrays.h"

#include "ah_param.h"

#include "ah_interpolator.h"
#include "ah_bicubic_interpolator.h"
#include "ah_sources.h"

void ah_integrator(const double H0, int *flag)
{
	// Previous values.
	double th_p, H_p, dH_p;

	// Intermediate values.
	double th, H, dH;

	// Weighted step.
	double dthw;

	// Integer counters.
	int k = 0, l = 0;

	// Interpolated values.
	double	i_a, i_b, i_h, i_c, i_phi,
		i_Aa, i_Ab, i_Ah, i_Ac, i_K,
		i_Dr_a, i_Dr_b, i_Dr_h, i_Dr_c, i_Dr_phi,
		i_Dz_a, i_Dz_b, i_Dz_h, i_Dz_c, i_Dz_phi;

	// RK4 sources arrays.
	double RKH[4] = {0.0, 0.0, 0.0, 0.0};
	double RKdH[4] = {0.0, 0.0, 0.0, 0.0};

	// Set initial values.
	array_H[0] = H0;
	array_dH[0] = 0.0;
	H = H0;
	dH = 0.0;
	th = 0.0;

	// Loop over NthTotal values.
	for (l = 0; l < NthTotal - 1; l++)
	{
		// Store old values.
		H_p = H;
		dH_p = dH;
		th_p = th;
		
		// Do RK4 integration.
		for (k = 1; k <= 4; k++)
		{
			// Find step size.
			if (k < 3)
			{
				dthw = 0.5 * dth;
			}
			else 
			{
				dthw = dth;
			}

			// Interpolate values.
			if (ahfind_bicubic)
			{
				ah_bicubic_interpolator(flag, H, th, 
					&i_a, &i_b, &i_h, &i_c, &i_phi,
					&i_Aa, &i_Ab, &i_Ah, &i_Ac, &i_K,
					&i_Dr_a, &i_Dr_b, &i_Dr_h, &i_Dr_c, &i_Dr_phi,
					&i_Dz_a, &i_Dz_b, &i_Dz_h, &i_Dz_c, &i_Dz_phi,
					a, Dr_a, Dz_a, Drz_a,
					b, Dr_b, Dz_b, Drz_b,
					h, Dr_h, Dz_h, Drz_h,
					c, Dr_c, Dz_c, Drz_c,
					phi, Dr_phi, Dz_phi, Drz_phi,
					A_a, Dr_A_a, Dz_A_a, Drz_A_a,
					A_b, Dr_A_b, Dz_A_b, Drz_A_b,
					A_h, Dr_A_h, Dz_A_h, Drz_A_h,
					A_c, Dr_A_c, Dz_A_c, Drz_A_c,
					K, Dr_K, Dz_K, Drz_K,
					Drr_a, Drrz_a,
					Drr_b, Drrz_b,
					Drr_h, Drrz_h,
					Drr_c, Drrz_c,
					Drr_phi, Drrz_phi,
					Dzz_a, Drzz_a,
					Dzz_b, Drzz_b,
					Dzz_h, Drzz_h,
					Dzz_c, Drzz_c,
					Dzz_phi, Drzz_phi);
			}
			else
			{
				ah_interpolator(flag, H, th, 
					&i_a, &i_b, &i_h, &i_c, &i_phi,
					&i_Aa, &i_Ab, &i_Ah, &i_Ac, &i_K,
					&i_Dr_a, &i_Dr_b, &i_Dr_h, &i_Dr_c, &i_Dr_phi,
					&i_Dz_a, &i_Dz_b, &i_Dz_h, &i_Dz_c, &i_Dz_phi,
					a,  b,  h,  c,  phi,
					A_a,  A_b,  A_h,  A_c,  K,
					Dr_a,  Dr_b,  Dr_h,  Dr_c, Dr_phi,
					Dz_a,  Dz_b,  Dz_h,  Dz_c, Dz_phi);
			}

			// Write interpolated values.
			if (k == 1)
			{
				array_a[l] = i_a;
				array_b[l] = i_b;
				array_h[l] = i_h;
				array_c[l] = i_c;
				array_phi[l] = i_phi;
			}

			// Check for error interpolation.
			if (*flag != 0)
			{
				return;
			}

			// Calculate sources.
			ah_sources(&RKH[k-1], &RKdH[k-1], th, H, dH,
				i_a, i_b, i_h, i_c, i_phi,
				i_Aa, i_Ab, i_Ah, i_Ac, i_K,
				i_Dr_a, i_Dr_b, i_Dr_h, i_Dr_c, i_Dr_phi,
				i_Dz_a, i_Dz_b, i_Dz_h, i_Dz_c, i_Dz_phi);
			
			// Update.
			if (k < 4)
			{
				H = H_p + dthw * RKH[k - 1];
				dH = dH_p + dthw * RKdH[k - 1];
				th = th_p + dthw;
			}
		}
		// RK4 integration.
		H = H_p + dth * (RKH[0] + 2.*RKH[1] + 2.*RKH[2] + RKH[3])/6.;
		dH = dH_p + dth * (RKdH[0] + 2.*RKdH[1] + 2.*RKdH[2] + RKdH[3])/6.;
		th = th_p + dth;

		// Write to memory.
		array_H[l + 1] = H;
		array_dH[l + 1] = dH;
	}

	// Obtain last values at equator.
	if (ahfind_bicubic)
	{
		ah_bicubic_interpolator(flag, H, th, 
			&i_a, &i_b, &i_h, &i_c, &i_phi,
			&i_Aa, &i_Ab, &i_Ah, &i_Ac, &i_K,
			&i_Dr_a, &i_Dr_b, &i_Dr_h, &i_Dr_c, &i_Dr_phi,
			&i_Dz_a, &i_Dz_b, &i_Dz_h, &i_Dz_c, &i_Dz_phi,
			a, Dr_a, Dz_a, Drz_a,
			b, Dr_b, Dz_b, Drz_b,
			h, Dr_h, Dz_h, Drz_h,
			c, Dr_c, Dz_c, Drz_c,
			phi, Dr_phi, Dz_phi, Drz_phi,
			A_a, Dr_A_a, Dz_A_a, Drz_A_a,
			A_b, Dr_A_b, Dz_A_b, Drz_A_b,
			A_h, Dr_A_h, Dz_A_h, Drz_A_h,
			A_c, Dr_A_c, Dz_A_c, Drz_A_c,
			K, Dr_K, Dz_K, Drz_K,
			Drr_a, Drrz_a,
			Drr_b, Drrz_b,
			Drr_h, Drrz_h,
			Drr_c, Drrz_c,
			Drr_phi, Drrz_phi,
			Dzz_a, Drzz_a,
			Dzz_b, Drzz_b,
			Dzz_h, Drzz_h,
			Dzz_c, Drzz_c,
			Dzz_phi, Drzz_phi);

	}
	else
	{
		ah_interpolator(flag, H, th, 
			&i_a, &i_b, &i_h, &i_c, &i_phi,
			&i_Aa, &i_Ab, &i_Ah, &i_Ac, &i_K,
			&i_Dr_a, &i_Dr_b, &i_Dr_h, &i_Dr_c, &i_Dr_phi,
			&i_Dz_a, &i_Dz_b, &i_Dz_h, &i_Dz_c, &i_Dz_phi,
			a,  b,  h,  c,  phi,
			A_a,  A_b,  A_h,  A_c,  K,
			Dr_a,  Dr_b,  Dr_h,  Dr_c, Dr_phi,
			Dz_a,  Dz_b,  Dz_h,  Dz_c, Dz_phi);
	}

	array_a[NthTotal - 1] = i_a;
	array_b[NthTotal - 1] = i_b;
	array_h[NthTotal - 1] = i_h;
	array_c[NthTotal - 1] = i_c;
	array_phi[NthTotal - 1] = i_phi;
}
