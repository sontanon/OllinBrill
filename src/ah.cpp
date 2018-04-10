#include "tools.h"

#include "param.h"
#include "arrays.h"

#include "ah_param.h"
#include "ah_tools.h"
#include "ah_geometry.h"
#include "ah_integrator.h"
#include "ah_mask.h"

void ah(void)
{
	// Auxiliary integer.
	int k = 0, l = 0;

	// Auxiliary double.
	double H0;

	// Auxiliary doubles that store the derivative value at the equator.
	double t1_dH, t2_dH;

	// Auxiliary integer used for integration success flag.
	int flag = 0;

	// Flag telling if last integration was succesfull.
	int last_int = 0;

	// Bisection method variables.
	double Ha, Hb, Hc;
	double dHa, dHb, dHc;

	// Exit tolerance.
	double tol = (ahfind_bicubic) ? pow(drr, 4) : pow(drr, 2);

	printf("APPARENT HORIZON FINDER IN AXISYMMETRY\n\n");

	printf(" ------------------ ------------------ \n");
	printf("|        H0        |     dH(pi/2)     |\n");
	printf("|------------------|------------------|\n");

	// Calculate extra geoemtry
	ah_geometry();

	// Initial values.
	k = 0;

	// Loop over values of rr.
	for (H0 = start + 20 * drr; H0 > ahfind_minr; H0 -= drr)
	{
		// First integrate.
		ah_integrator(H0, &flag);

		// Check for integration success.
		if (flag != 0)
		{
			// Print that integration exceeded grid extension.
			printf("|    %-6.2E      |      OUTSIDE     |\n", H0);

			// Reset flag to 0.
			flag = 0;
			
			// Indicate that last integration failed.
			last_int = 0;

			// Cycle.
			continue;
		}

		// If integration succeeded we store the value of dH(pi/2) at
		// temporary variable t1_dH.
		t1_dH = array_dH[NthTotal - 1];

		// Print to screen.
		printf("|    %-6.2E      |    %-6.2E      |\n", H0, t1_dH);

		// Check if we can compare to last integration.
		if (last_int && (t1_dH * t2_dH <= 0.0))
		{
			// Check if we have sign change. If so, we have found an
			// apparent horizon.
			// Print to screen that AH was found.
			printf(" ------------------ ------------------ \n");
			printf("\nAPPARENT HORIZON WAS FOUND!\n\n");

			// Set found flag to 1.
			ah_found = 1;

			// Now comes the bisection refinement.
			Ha = H0;
			Hb = H0 + drr;
			dHa = t1_dH;
			dHb = t2_dH;

			for (l = 0; l < ahfind_bis_iter; l++)
			{
				// Get midpoint.
				Hc = 0.5 * (Ha + Hb);

				// Integrate.
				ah_integrator(Hc, &flag);

				// Check for integration success.
				if (flag != 0)
				{
					printf("\nBISSECTION REFINEMENT FAILED AT INTEGRATION!\n");
					return;
				}

				// Get derivative at equator.
				dHc = array_dH[NthTotal - 1];

				// Check for convergence.
				if (ABS(Hc - Ha) < tol)//|| ABS(dHc) < tol)
				{
					// Exit.
					printf("\nAPPARENT HORIZON WAS REFINED USING %d BISSECTION ITERATIONS!\n\n", l + 1);
					// Write data to files.
					write_ah(t);

					// Find maximum and minimum.
					ah_max = fabs(array_H[cblas_idamax(NthTotal, array_H, 1)]);
					ah_min = fabs(array_H[cblas_idamin(NthTotal, array_H, 1)]);

					// Calculate mask.
					ah_mask();

					start = Hc;

					return;
				}

				// Examine sign of three points to cycle.
				if (dHc * dHa <= 0)
				{
					// Zero is in [a,c].
					Hb = Hc;
					dHb = dHc;
				}
				else
				{
					// Zero is in [c,b].
					Ha = Hc;
					dHa = dHc;
				}
			}
			printf("\nBISSECTION REFINEMENT FAILED AFTER %d MAX ITERATIONS.\n", ahfind_bis_iter);

			// Write data to files.
			write_ah(t);

			// Find maximum and minimum.
			ah_max = fabs(array_H[cblas_idamax(NthTotal, array_H, 1)]);
			ah_min = fabs(array_H[cblas_idamin(NthTotal, array_H, 1)]);

			// Calculate mask.
			ah_mask();

			return;
		}
		// If the last integration was unsuccessfull, we cannot compare.
		else
		{
			// Set last_int to 1 because this integration was
			// succesfull.
			last_int = 1;

			// Cycle variables.
			t2_dH = t1_dH;
		}
	}

	// If we reach this point, then there we found no apparent horizon in our data.
	printf(" ------------------ ------------------ \n");
	printf("\nNO APPARENT HORIZON WAS FOUND IN THE INTEGRATION RANGE.\n\n");
	return;
}
