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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <gmp.h>
//#include <omp.h>

#ifdef WITH_PARI
#include <pari/pari.h>
#endif

#include "verify.h"
#include "functions.h"
#include "sieve.h"

int main(int argc, char ** argv)
{
	if (argc != 4)
	{
		printf("Format: ./verify [D_max] [files] [folder]\n");
		exit(1);
	}

	long D_max = atol(argv[1]);
	long files = atol(argv[2]);
	const char * folder = argv[3];

	if (D_max < 0)
	{
		perror("D_max should should be non-negative.\n");
		exit(1);
	}

	long D_sqrt = sqrt(D_max) + 1;

	int * primes = (int *) malloc(((unsigned int) (1.25506 * D_sqrt / log(D_sqrt))) * sizeof(int));
	prime_sieve(D_sqrt, primes);
	primes = (int *) realloc(primes, (2 + primes[0]) * sizeof(int));

	printf("\n\n%ld discriminants, %ld files, %d threads\n", D_max, files, /*omp_get_max_threads()*/1);
	fflush(stdout);

	#ifdef WITH_PARI
	pari_init(100000, 0);
	#endif

	mpz_t LHS;
	mpz_init(LHS);
	left_hand_side(LHS, D_max, files, primes, folder);

	gmp_printf("Left hand side computed.\n%Zd\n\n", LHS);
	fflush(stdout);

	#ifdef WITH_PARI
	pari_close();
	#endif

	mpz_t RHS;
	mpz_init(RHS);
	right_hand_side(RHS, D_max / 8, MIN(D_max / 8, FAC_TOTAL), primes);

	gmp_printf("Right hand side computed.\n%Zd\n\n", RHS);
	fflush(stdout);

	int result = mpz_cmp(LHS, RHS);

	if (result == 0)
	{
		printf("Sides are equal\n");
	}
	else if (result > 0)
	{
		printf("Left operand is bigger.\n");
	}
	else if (result < 0)
	{
		printf("Right operand is bigger\n");
	}

	mpz_clear(LHS);
	mpz_clear(RHS);

	return 0;
}
