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

void partial_left_hand_side(mpz_t LHS, const char * file, const char * folder, const int index, const long blocksize, const long x, const int a, const int m, const int * primes, int ** factors)
{
	char name[500];
	
	sprintf(name, "gunzip -c %s/%s%d.gz > %s/file%d", folder, file, index, folder, index);
	system(name);

	char kron[primes[0] + 1];

#ifdef DEBUG
printf("%d mod %d, n = %ld\n", a, m, x);
#endif

	sprintf(name, "%s/file%d", folder, index);
	FILE * clgrp = fopen(name, "r");
	
	char line[500];
	char * end;
	
	long dist, h, h_omega, r, K, Y, bound, i, j, k, l, t, total, rem;
	long D = index * blocksize + a;
	long p, D_next;
	int * i_factors, * t_factors;
	long div[1000];

	mpz_t sum;
	mpz_init(sum);

	while (fgets(line, 500, clgrp) != NULL)
	{
		end = NULL;
		dist = strtol(line, &end, 10);
		D += dist * m;
		h = strtol(end, &end, 10);

#ifdef DEBUG
printf("h(-%ld) = %ld, ", D, h);
#endif

		// compute H(Delta)
		h_omega = (D > 4) ? (6 * h) : ((D == 4) ? 3 : 2);

		rem = (D & 7);
		if (rem != 3)
		{
			Y = ((long) sqrt((x << 3) - D));
			if (rem == 7)
			{
				r = (Y + 1) / 2;
			}
			else
			{
				Y &= -2;	// last even
				r = (Y >> 2) + ((Y & 3) == 2 && rem == 4);
			}

			mpz_add_ui(sum, sum, 2 * r * h_omega);

#ifdef DEBUG
printf("H(-%ld) = %ld, r = %ld\n", D, 2 * r * h_omega * K, r);
#endif

			#ifdef WITH_PARI
			GEN g = stoi(D);
			GEN H = hclassno(g);

			long H_pari = (typ(H) == t_INT) ? (6 * gtolong(H)) : (6 * gtolong(gel(H, 1)) / gtolong(gel(H, 2)));

			if (H_pari != h_omega * K)
			{
				printf("PARI ERROR FOR D=%ld: H_pari=%ld, H=%ld\n", D, H_pari, h_omega * K);
				exit(1);
			}

			cgiv(H);
			cgiv(g);
			#endif

			// account for H(-8n)
			if (rem == 0)
			{
#ifdef DEBUG
printf("A Account for H(-%ld): %ld, h_omega=%ld, K=%ld\n", D, h_omega * K, h_omega, K);
#endif
				mpz_add_ui(sum, sum, h_omega);
			}
		}
		
		bound = sqrt((x << 3) / ((double) D));
		
		for (i = 2, j = 4; i <= bound; j += (i << 1) + 1, i++)
		{
			// compute new discriminant
			D_next = D * j;
			
			rem = (D_next & 7);
			if (rem != 3)
			{
				Y = ((long) sqrt((x << 3) - D_next));
				if (rem == 7)
				{
					r = (Y + 1) >> 1;
				}
				else
				{
					Y &= -2;	// last even
					r = (Y >> 2) + (rem && (Y & 3) == 2);
				}

				// compute all kronecker symbols
				i_factors = factors[i];

				for (k = 1; k <= i_factors[0]; k++)
				{
					kron[i_factors[k]] = kronecker_symbol(-D, primes[i_factors[k]]);
				}

				// compute K
				total = divisors_list(div, i, i_factors, primes, WITH_INDICES);

				K = 1;
				for (k = 1; k < total; k++)
				{
					t = div[k];
					t_factors = factors[t];
#ifdef DEBUG
printf("%ld, #factors=%d\n", t, t_factors[0]);
#endif

					for (l = 1; l <= t_factors[0]; l++)
					{
						p = primes[t_factors[l]];
#ifdef DEBUG
printf("(%ld,%d)=%d\n", -D, p, kron[t_factors[l]]);
#endif
						t /= p;
						t *= (p - kron[t_factors[l]]);
					}

					K += t;
				}

				mpz_add_ui(sum, sum, 2 * r * h_omega * K);

#ifdef DEBUG
printf("H(-%ld=%ld*%ld^2) = %ld, K=%ld, Y=%ld, r = %ld\n", D_next, D, i, r * 2 * h_omega * K, K, Y, r);
#endif

				#ifdef WITH_PARI
				GEN g = stoi(D_next);
				GEN H = hclassno(g);

				long H_pari = (typ(H) == t_INT) ? (6 * gtolong(H)) : (6 * gtolong(gel(H, 1)) / gtolong(gel(H, 2)));

				if (H_pari != h_omega * K)
				{
					printf("PARI ERROR FOR D=%ld: H_pari=%ld, H=%ld\n", D_next, H_pari, h_omega * K);
					exit(1);
				}

				cgiv(H);
				cgiv(g);
				#endif

				// account for H(-8n)
				if (rem == 0)
				{
#ifdef DEBUG
printf("B Account for H(-%ld): %ld, h_omega=%ld, K=%ld\n", D_next,  h_omega * K, h_omega, K);
#endif
					mpz_add_ui(sum, sum, h_omega * K);
				}
			}
		}
	}

	mpz_add(LHS, LHS, sum);
	mpz_clear(sum);

	sprintf(name, "%s/file%d", folder, index);	
	fclose(clgrp);
	remove(name);
}

void left_hand_side(mpz_t LHS, const long D_max, const long files, const int * primes, const char * folder)
{
	struct timeval begin, end;
	unsigned long exec_time;

	mpz_set_ui(LHS, 0);

	long size = 1;
	long temp = 1;
	long D_sqrt = sqrt(D_max) + 1;

	while (temp < D_sqrt)
	{
		temp *= primes[size++];
		size++;
	}

	int ** factors = (int **) malloc(D_sqrt * sizeof(int *));

	for (long i = 0; i < D_sqrt; i++)
	{
		factors[i] = (int *) malloc(size * sizeof(int));
	}

	regular_sieve(D_sqrt, D_sqrt, factors, primes, WITH_INDICES);

	int a[4] = {4, 8, 3, 7};
	int m[4] = {16, 16, 8, 8};
	char * file[4] = { "cl4mod16/cl4mod16.", "cl8mod16/cl8mod16.", "cl3mod8/cl3mod8.", "cl7mod8/cl7mod8." };

	mpz_t mod_sum;
	mpz_init(mod_sum);

	for (int i = 0; i < 4; i++)
	{
		gettimeofday(&begin, NULL);
		mpz_set_ui(mod_sum, 0);

		#pragma omp parallel for
		for (int j = 0; j < files; j++)
		{
			mpz_t sum;
			mpz_init(sum);
			partial_left_hand_side(sum, file[i], folder, j, D_max / files, D_max / 8, a[i], m[i], primes, factors);
			
			#pragma omp critical
			{
				mpz_add(mod_sum, mod_sum, sum);
			}
			mpz_clear(sum);
		}

#ifdef DEBUG
gmp_printf("Sum %d mod %d: %Zd\n\n", a[i], m[i], mod_sum);
#endif

		mpz_add(LHS, LHS, mod_sum);
		gettimeofday(&end, NULL);
		exec_time = (end.tv_sec * 1e6 + end.tv_usec) - (begin.tv_sec * 1e6 + begin.tv_usec);
		printf("%d mod %d processed in %.3f\n", a[i], m[i], exec_time / 1e6);
		gmp_printf("%Zd\n", mod_sum);
		fflush(stdout);
	}

	mpz_clear(mod_sum);

#ifdef DEBUG
gmp_printf("LHS = %Zd\n\n", LHS);
#endif

	for (long i = 0; i < D_sqrt; i++)
	{
		free(factors[i]);
	}

	free(factors);
}

void partial_right_hand_side(mpz_t sum, const long blocksize, const long l, const int * primes, int ** factors)
{
	segmented_sieve(sqrt((l + blocksize) << 1), blocksize, l, factors, primes, 0);

	int * f;

	if (l == 0)
	{
		factors[1][1] = 2;
	}

	for (long i = (l == 0) ? 2 : ((l & 1) ? 0 : 1); i < blocksize; i += 2)
	{
		f = factors[i];
		f[++f[0]] = 2;
	}

	const long  max = (l + blocksize) << 1;

	long i, n, root, total;
	long div[MAX_DIVISORS];

	unsigned long psi;

	char is_square;

	long cur_f = (l == 0) ? 1 : 0;

	for (n = (l == 0) ? 2 : (l << 1); n < max; cur_f++, n += 2)
	{
		psi = n;
		total = divisors_list(div, n, factors[cur_f], primes, 0);

		root = sqrt(n);
		is_square = (root * root == n);

		for (i = 1; i < total - 1; i++)
		{
			if (div[i] > root)
			{
				psi += div[i];
			}
		}

		psi *= 12;

		if (is_square)
		{
			psi += 6 * root + 1;
		}

#ifdef DEBUG
printf("%ld: %lu\n", n, psi);
#endif

		mpz_add_ui(sum, sum, psi);
	}
}

void right_hand_side(mpz_t RHS, const long n_max, const long blocksize, const int * primes)
{
	struct timeval begin, end;
	unsigned long exec_time;

	gettimeofday(&begin, NULL);
	mpz_set_ui(RHS, 0);

	long temp = 1;
	long size = 1;

	while (temp < (n_max << 1))
	{
		temp *= primes[size++];
		size++;
	}

	const int max_threads = /*omp_get_max_threads()*/1;

	int ** factors[max_threads];

	for (int n = 0; n < max_threads; n++)
	{
		factors[n] = (int **) malloc(blocksize * sizeof(int *));

		for (long i = 0; i < blocksize; i++)
		{
			factors[n][i] = (int *) malloc((size + 1) * sizeof(int));
		}
	}

	const long num_blocks = n_max / blocksize;

#ifdef DEBUG
printf("# blocks: %d\n", num_blocks);
#endif

	#pragma omp parallel for
	for (int i = 0; i < num_blocks; i++)
	{
		mpz_t sum;
		mpz_init(sum);
		partial_right_hand_side(sum, blocksize, i * blocksize + 1, primes, factors[/*omp_get_thread_num()*/0]);
		
		#pragma omp critical
		{
			mpz_add(RHS, RHS, sum);
		}
		mpz_clear(sum);
	}

#ifdef DEBUG
gmp_printf("RHS = %Zd\n\n", RHS);
#endif

	for (int n = 0; n < max_threads; n++)
	{
		for (long i = 0; i < blocksize; i++)
		{
			free(factors[n][i]);
		}

		free(factors[n]);
	}
	gettimeofday(&end, NULL);
	exec_time = (end.tv_sec * 1e6 + end.tv_usec) - (begin.tv_sec * 1e6 + begin.tv_usec);
	printf("Right hand side processed in %.3f\n", exec_time / 1e6);
}
