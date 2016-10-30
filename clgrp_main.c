/*=============================================================================

    This file is part of CLGRP.

    CLGRP is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    CLGRP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CLGRP; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2016 Anton Mosunov
 
******************************************************************************/

#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include "clgrp.h"
#include "sieve.h"

int main(int argc, char * argv[])
{
	MPI_Init(&argc, &argv);

	if (argc != 7)
	{
		perror("Format: mpirun -np [#procs] ./clgrp [D_max] [files] [a] [m] [h_prefix] [folder] [D_total]\nSet h_prefix to \"null\" if class numbers were not precomputed.\n");
		exit(1);
	}

	long D_max = atol(argv[1]);
	const long files = atol(argv[2]);
	const int a = atoi(argv[3]);
	const int m = atoi(argv[4]);
	const char * h_prefix = argv[5];
	const char * folder = argv[6];
	const long D_total = D_max / (files * m);

	int * primes;
	int ** h_factors;

	if (strcmp(h_prefix, "null") == 0)
	{
		h_prefix = NULL;
	}

	int i, myrank, idx = 0;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if (myrank == 0)
	{
		int num_procs;
		MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
		num_procs--;

		for (i = 0; i < num_procs; i++)
		{
			MPI_Send(&i, 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD);
		}

		for (i = num_procs; i < files; i++)
		{
			MPI_Recv (&idx, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(&i, 1, MPI_INT, idx, 0, MPI_COMM_WORLD);
		}

		for (i = 0, idx = -1; i < num_procs; i++)
		{
			MPI_Send(&idx, 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD);
		}
	}
	else
	{
		#ifdef WITH_PARI
		pari_init(1000000, 0);
		#endif

		long D_root = sqrt(D_max);

		primes = (int *) malloc(((int) (1.25506 * D_root / log(D_root))) * sizeof(int));
		prime_sieve(D_root, primes);
		primes = (int *) realloc(primes, (2 + primes[0]) * sizeof(int));

		int h_max = 0;
		long temp = 1;

		int * h_list = NULL;

		if (h_prefix)
		{
			h_list = (int *) malloc(FAC_TOTAL * sizeof(int));

			// Ramare's bound
			h_max = (1/M_PI) * sqrt(D_max) * (0.5 * log(D_max) + 2.5 - log(6)) + 1;

			// compute maximal number of prime factors of a class number
			for (i = 1; temp < h_max; i++)
			{
				temp *= primes[i];
			}

			const int h_max_factors = i;
			h_factors = (int **) malloc(h_max * sizeof(int *));

			for (i = 0; i < h_max; i++)
			{
				h_factors[i] = (int *) malloc(h_max_factors * sizeof(int));
			}

			regular_sieve(h_max, h_max, h_factors, primes, 0);
		}

		// compute maximal number of prime factors of a discriminant
		temp = ((a & 3) == 0) ? 4 : 1;

		for (i = 2; temp < D_max; i++)
		{
			temp *= primes[i];
		}

		const int size = i - ((a & 3) != 0) + 1;

		int ** factors = (int **) malloc(FAC_TOTAL * sizeof(int *));

		for (i = 0; i < FAC_TOTAL; i++)
		{
			factors[i] = (int *) malloc(size * sizeof(int));
		}

		MPI_Recv(&idx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		while (idx != -1)
		{
			tabulate_bjt(idx, D_total, h_prefix, folder, a, m, primes, h_factors, factors, h_list);
			MPI_Send(&myrank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Recv(&idx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		free(primes);

		if (h_prefix)
		{
			for (long i = 0; i < h_max; i++)
			{
				free(h_factors[i]);
			}

			free(h_factors);
			free(h_list);
		}


		for (int i = 0; i < FAC_TOTAL; i++)
		{
			free(factors[i]);
		}

		free(factors);

		#ifdef WITH_PARI
		pari_close();
		#endif
	}

	MPI_Finalize();

	return 0;
}
