// Headers.
#include "tools.h"
#include "parser.h"
#include "allocate_arrays.h"
#include "array_init.h"
#include "evolve.h"
#include "write_files.h"

// Global variables.
#define MAIN_FILE
#include "param.h"
#include "arrays.h"

// FILES.
#define FILE_MAIN_FILE
#include "files.h"

// WRITE FLAGS.
#define FLAGS_MAIN_FILE
#include "write_flags.h"

// DISPLAY NUMBER OF OMP THREADS.
#define OMP_DEBUG

int main(int argc, char *argv[])
{
	// Auxiliary integers.
	int i, j, k;

	// Auxiliary doubles.
	double aux_r, aux_z, aux_rr;

	// Initial message.
	printf("**************************************\n");
	printf("**************************************\n");
	printf("***                                ***\n");
	printf("***           OLLINBRILL           ***\n");
	printf("***                                ***\n");
	printf("***                                ***\n");
	printf("***     First Revision: 6/2017     ***\n");
	printf("***                                ***\n");
	printf("***     Last  Revision: 3/2018     ***\n");
	printf("***                                ***\n");
	printf("**************************************\n");

	// File name is in argv[1]. Check that we have at
	// least one argument.
	if (argc < 2)
	{
		printf("***                                ***\n");
		printf("***  Usage: ./OllinBrill file.par  ***\n");
		printf("***                                ***\n");
		printf("***    Missing parameter file.     ***\n");
		printf("***                                ***\n");
		printf("**************************************\n");
		printf("**************************************\n");
		return(EXIT_FAILURE);
	}

	// Parse arguments.
	parser(argv[1]);

	// FOR NOW THE CODE IS LIMITED TO dr == dz.
	assert(dr == dz);

	// Create output directory.
	char opt;
	struct stat st = { 0 };
	if (stat(dirname, &st) == -1)
	{
#ifdef WIN
		_mkdir(dirname);
#else
		mkdir(dirname, 0755);
#endif
	}
	else
	{
		printf("***                                ***\n");
		printf("***   WARNING: Output directory    ***\n");
		printf("***            %-15s    ***\n", dirname);
		printf("***            already exists!     ***\n");
		printf("***                                ***\n");
		printf("***   Press (y/n) to proceed: ");
		opt = getchar();
		printf("***                                ***\n");
		if ((opt == 'y') || (opt == 'Y'))
		{
			printf("***     User chose to continue.    ***\n");
			printf("***                                ***\n");
		}
		else
		{
			printf("***       User chose to abort.     ***\n");
			printf("***                                ***\n");
			printf("**************************************\n");
			printf("**************************************\n");
			return(EXIT_FAILURE);
		}
	}


	// Copy parameter file to directory.
	char cmd[256];
	memset(cmd, 0, 256);
#ifdef WIN
	sprintf(cmd, "COPY %s %s", argv[1], dirname);
#else
	sprintf(cmd, "cp %s %s", argv[1], dirname);
#endif
	system(cmd);

	// Cd to output directory.
	if (chdir(dirname) == -1)
	{
		printf("***                                ***\n");
		printf("***   WARNING: Could not cd to     ***\n");
		printf("***            output directory!   ***\n");
		printf("***                                ***\n");
		printf("***    Press (y/n) to write in     ***\n");
		printf("***       current directory: ");

		opt = getchar();
		printf("***                                ***\n");
		if ((opt == 'y') || (opt == 'Y'))
		{
			printf("***     User chose to continue.    ***\n");
			printf("***                                ***\n");
		}
		else
		{
			printf("***       User chose to abort.     ***\n");
			printf("***                                ***\n");
			printf("**************************************\n");
			printf("**************************************\n");
			return(EXIT_FAILURE);
		}
	}

	// Print main arguements.
	printf("**************************************\n");
	printf("***                                ***\n");
	printf("***   Evolving Brill waves using   ***\n");
	printf("***    BSSN Numerical Relativity   ***\n");
	printf("***                                ***\n");
	printf("***   dirname    = %-15s ***\n", dirname);
	printf("***                                ***\n");
	printf("***   Slicing    = %-13s   ***\n", slicing);
	printf("***   ilapse     = %-13s   ***\n", ilapse);
	printf("***   Shift      = %-13s   ***\n", shift);
	if (!(strcmp(shift, "none") == 0))
		printf("***   ishift     = %-13s   ***\n", ishift);
	printf("***   Flavor     = %-13s   ***\n", bssn_flavor);
	printf("***   Integrator = %-13s   ***\n", integrator);
	printf("***   Order      = %-13s   ***\n", order);
	printf("***   Boundary   = %-13s   ***\n", boundtype);
	printf("***                                ***\n");
	printf("***     Space step = %-5.3E     ***\n", dr);
	printf("***     Final time = %-5.3E     ***\n", FinalTime);
	printf("***                                ***\n");
	printf("**************************************\n");

#ifdef OMP_DEBUG
#pragma omp parallel
	{
#pragma omp master
		{
			// Determine OMP threads.
			printf("**************************************\n");
			printf("***                                ***\n");
			printf("***    Maximum OMP threads = %d     ***\n", omp_get_max_threads());
			printf("***    Currently running on %d      ***\n", omp_get_num_threads());
			printf("***                                ***\n");
			printf("**************************************\n");
		}
	}
#endif

	// Determine sigma via BSSN flavor (irrelevant for no shift).
	if (strcmp(bssn_flavor, "lagrangian") == 0)
	{
		sigma = 1.0;
	}
	else if (strcmp(bssn_flavor, "eulerian") == 0)
	{
		sigma = 0.0;
		printf("\n\nBSSN EULERIAN FLAVOR NOT SUPPORTED YET!\n");
		return(EXIT_FAILURE);
	}

	// Determine ghost zones according to code order.
	if (strcmp(order, "two") == 0)
	{
		ghost = 2;
	}
	else if (strcmp(order, "four") == 0)
	{
		ghost = 3;
	}

	printf("**************************************\n");
	printf("***                                ***\n");
	printf("***         Ghost zones = %1d        ***\n", ghost);
	printf("***                                ***\n");
	printf("**************************************\n");

	// Determine linear dimension.
	NrTotal = NrInterior + ghost + 1;
	NzTotal = NzInterior + ghost + 1;
	DIM = NrTotal * NzTotal;

	// Allocate arrays.
	allocate_arrays(true);

	// Set all memory to zero.
	array_init();

	printf("**************************************\n");
	printf("***                                ***\n");
	printf("***   Allocated arrays of total    ***\n");
	printf("***   number of points = %-9d ***\n", DIM);
	printf("***                                ***\n");
	printf("**************************************\n");

	// Fill r, z and rr grids.
#pragma omp parallel shared(r, z, rr, r2) private(aux_r, aux_z, j)
	{
#pragma omp for schedule(guided)
		// Loop in r direction.
		for (i = 0; i < NrTotal; i++)
		{
			aux_r = ((double)(i - ghost) + 0.5) * dr;

			// Loop in z direction.
			for (j = 0; j < NzTotal; j++)
			{
				aux_z = ((double)(j - ghost) + 0.5) * dz;
				r[IDX(i, j)] = aux_r;
				r2[IDX(i, j)] = aux_r * aux_r;
				z[IDX(i, j)] = aux_z;
				rr[IDX(i, j)] = sqrt(aux_r * aux_r + aux_z * aux_z);
			}
		}
	}

	// Write coordinate grids.
	writeSingleFile(r, "r.asc");
	writeSingleFile(z, "z.asc");

	// Open files.
	open_files();

	// Evolution.
	evolve();

	// Deallocate arrays.
	allocate_arrays(false);

	// Close files.
	close_files();

	// Destroy configuration.
	config_destroy(&cfg);

	printf("**************************************\n");
	printf("***                                ***\n");
	printf("***       Deallocated arrays       ***\n");
	printf("***                                ***\n");
	printf("**************************************\n");

	// Print final message.
	printf("**************************************\n");
	printf("***                                ***\n");
	printf("***   All done! Have a nice day!   ***\n");
	printf("***                                ***\n");
	printf("**************************************\n");
	printf("**************************************\n");

	return 0;
}