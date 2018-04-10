#include "tools.h"

#include "param.h"
#include "arrays.h"

#include "sources_geometry.h"
#include "sources_lapse.h"
#include "shift_sources.h"
#include "accumulate.h"
#include "boundary.h"
#include "update.h"
#include "symmetries.h"
#include "auxiliary_geometry.h"
#include "gauge.h"
#include "shift_geometry.h"
#include "write_files.h"
#include "constraints.h"

void one_step(const double dt)
{
	// Extra variables.
	int k, niter;
	double dtw, weight;
	const double third = 1. / 3.;
	const double sixth = 1. / 6.;

	// FIND NUMBER OF INTERNAL ITERATIONS.
	if (strcmp(integrator, "icn") == 0)
	{
		niter = 3;
	}
	else if (strcmp(integrator, "rk4") == 0)
	{
		niter = 4;
	}

	// START INTERNAL ITERATIONS.
	for (k = 1; k <= niter; k++)
	{
		// FIND WEIGHTS.
		// ICN.
		if (strcmp(integrator, "icn") == 0)
		{
			if (k < 3)
			{
				dtw = 0.5 * dt;
			}
			else
			{
				dtw = dt;
			}
		}
		// RK4.
		else if (strcmp(integrator, "rk4") == 0)
		{
			switch (k)
			{
			case 1:
				dtw = 0.5 * dt;
				weight = sixth;
				break;
			case 2:
				dtw = 0.5 * dt;
				weight = third;
				break;
			case 3:
				dtw = dt;
				weight = third;
				break;
			case 4:
				dtw = dt;
				weight = sixth;
				break;
			}
		}

		// Calculate geometry sources.
		sources_geometry();

		// If not maximal slicing, calculate sources for the lapse.
		sources_lapse();

#ifdef SHIFT
		// If shift, calculate sources for shift.
		shift_sources();
#endif

		// Boundary conditions.
		boundary();

		// RUNGE-KUTTA ACCUMULATE.
		if (strcmp(integrator, "rk4") == 0)
		{
			accumulate(k, niter, weight);
		}

		// UPDATE VARIABLES.
		update(dtw);

		// SYMMETRIES.
		symmetries();

		// SPACETIME
		// Update geometry.
		auxiliary_geometry();

		// Update gauge.
		if (strcmp(slicing, "maximal") == 0)
		{
			// Maximal slicing does not have to be solved at 
			// every iteration.
			if (strcmp(integrator, "icn") == 0)
			{
				// If we have ICN integration, we only have to
				// calculate the gauge at t = t0 + 0.5 * dt and
				// t = t0 + dt which correspond to iterations
				// k = 1 and k = 3. k = 2 corresponds to 
				// t = t0 + 0.5 * dt which was calculated on k = 1.
				if (strcmp(ell_freq, "every") == 0)
				{
					gauge();
				}
				else if (strcmp(ell_freq, "half") == 0)
				{
					if (k == 1 || k == 3)
					{
						gauge();
					}
				}
				else if (strcmp(ell_freq, "full") == 0)
				{
					if (k == 3)
					{
						gauge();
					}
				}
			}
			else if (strcmp(integrator, "rk4") == 0)
			{
				// In a similar manner, RK4 only needs lapse 
				// calculation at steps 1 and 4.
				if (strcmp(ell_freq, "every") == 0)
				{
					gauge();
				}
				else if (strcmp(ell_freq, "half") == 0)
				{
					if (k == 1 || k == 4)
					{
						gauge();
					}
				}
				else if (strcmp(ell_freq, "full") == 0)
				{
					if (k == 4)
					{
						gauge();
					}
				}
			}
		}
		// Calculate gauge for non-maximal slicing.
		else
		{
			gauge();
		}

#ifdef SHIFT
		// Shift terms.
		if (!(strcmp(shift, "none") == 0))
		{
			shift_geometry();
		}
#endif
	}

	return;
}
