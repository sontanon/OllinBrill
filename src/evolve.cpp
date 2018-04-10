#include "tools.h"

#include "param.h"
#include "arrays.h"

#include "initial.h"
#include "symmetries.h"

#include "auxiliary_geometry.h"
#include "gauge.h"
#include "shift_geometry.h"
#include "analysis_metric.h"

#include "write_files.h"
#include "constraints.h"

#include "save_old.h"
#include "one_step.h"

#include "pardiso_init.h"
#include "pardiso_die.h"
#include "low_rank.h"

#ifdef AH
#include "ah_tools.h"
#include "ah.h"
#endif

#undef DEBUG_ONE

void evolve(void)
{
	// Full dr step counter.
	int l = 0;

	double zero = 0.0;

	// Initialize time.
	t = zero;

	// Find initial data.
	initial();

	// Assert r, z symmetries.
	symmetries();

	// Update auxiliary geometry.
	auxiliary_geometry();

	// Update gauge functions.
	gauge();

#ifdef SHIFT
	// Calculate shift terms.
	if (!(strcmp(shift, "none") == 0))
	{
		shift_geometry();
	}
#endif

	// Do metric analysis, i.e. CFL condition.
	metric_analysis();

	// Set initial Courant factor and time step.
	dtfac_old = dtfac_new;
	CFL_n_old = CFL_n_new;
	CFL_count = 0;
	CFL_stop = 2 << (CFL_n_old - 1);
	double dt = dtfac_old * dr;

	printf("***                                ***\n");
	printf("***         dtfac = %-5lf       ***\n", dtfac_old);
	printf("***                                ***\n");


#ifdef AH
	// AH finder.
	if (ahfind_flag)
	{
		// First initialize.
		ah_init();

		printf("***                                ***\n");
		printf("***     Initialized AH finder      ***\n");
		printf("***                                ***\n");
		printf("**************************************\n");

		// Next, look for AH in initial data.
		ah();
	}
#endif

	// Initial constraints.
	constraints(0);

	// Save initial data.
	write_files(t);

	// Output to screen.
	printf("**************************************\n");
	printf("***                                ***\n");
	printf("***   Saved initial data to files  ***\n");
	printf("***                                ***\n");
	printf("**************************************\n");
	printf("**************************************\n");
	printf("***                                ***\n");
	printf("***   Starting main evolution...   ***\n");
	printf("***                                ***\n");

	// Initialize maximal slicing solver.
	if (strcmp(slicing, "maximal") == 0)
	{
		// Initialize memory.
		pardiso_init(NrInterior, NzInterior);

		// Set permutation and flag and counter.
		perm_flag = (use_permutation) ? 2 : 0;
		perm_count = 0;

		// Preconditioner flag is set to 0 whether or not we utilize it.
		precond_flag = 0;

		// Low rank initialization and setup.
		if (use_low_rank)
		{
			// No permutation or preconditioner may be used.
			use_permutation = 0;
			use_preconditioner = 0;
			
			// Auxiliary numerical order.
			int norder;
			if (strcmp(order, "two") == 0)
				norder = 2;
			else if (strcmp(order, "four") == 0)
				norder = 4;

			low_rank_init(NrInterior, NzInterior, norder);

			low_rank_set(NrInterior, NzInterior, norder);

			// Set initial flag to zero.
			low_rank_flag = 0;
		}
	}

	// Start main evolution loop.
	while (t < FinalTime)
	{
		// EVOLVE.
		// Time.
		t += dt;

		// CFL count.
		CFL_count += 1;

		printf("***                                ***\n");
		printf("***           t = %-5lf         ***\n", t);
		printf("***                                ***\n");

		// Save old data.
		save_old();

		// Advance one time step via RK4 or ICN integrators.
		one_step(dt);

		// EVOLVE IS DONE.
		// Calculate dtfac: now stored in new variables.
		metric_analysis();

		// Check if counter has done full dr step.
		if (CFL_count == CFL_stop)
		{
			// Increment l counter.
			l += 1;

			// Reset dtfac.
			dtfac_old = dtfac_new;
			CFL_n_old = CFL_n_new;
			CFL_count = 0;
			CFL_stop = 2 << (CFL_n_old - 1);
			dt = dtfac_old * dr;

#ifdef AH
			// Check for AH conditions.
			if (ahfind_flag && ((l % ahfind_freq) == 0) && (t >= ahfind_start_time))
			{
				ah();
			}
#endif
			// Constraint analysis every constraintCheck steps.
			if ((l % constraintCheck) == 0)
			{
				// Calculate constraints.
				constraints(l);
			}
			// Save data.
			if ((l % Noutput2D) == 0)
			{
				write_files(t);
			}
		}
		// Check if new CFL factor is bigger than previous. If so, we remain with previous.
		if (dtfac_old <= dtfac_new)
		{
			dtfac_new = dtfac_old;
			CFL_n_new = CFL_n_old;
		}
		// Check if new CFL factor is smaller than previous. If so, we have to refine time step.
		else
		{
			dtfac_old = dtfac_new;
			CFL_count *= 2 << (CFL_n_new - CFL_n_old - 1);
			CFL_n_old = CFL_n_new;
			CFL_stop = 2 << (CFL_n_new - 1);
			dt = dtfac_old * dr;
		}
	}

	// Deallocate maximal slicing solver.
	pardiso_die();

	// Deallocate low rank structures.
	if (use_low_rank)
	{
		low_rank_die();
	}

#ifdef AH
	// Deallocate AH structures.
	if (ahfind_flag)
	{
		ah_destroy();
	}
#endif

	printf("***                                ***\n");
	printf("***  Evolution loop has finished!  ***\n");
	printf("***                                ***\n");
	printf("**************************************\n");

	return;
}
