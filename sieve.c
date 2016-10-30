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
#include <math.h>
#include <dirent.h>

#include <liboptarith/gcd/gcd_binary_l2r.h>

#include "functions.h"
#include "sieve.h"


void prime_sieve(const int max_prime, int * primes)
{
	long i, j;

	primes[0] = 0;

	char is_prime[max_prime];
	memset(is_prime, 1, max_prime);

	for (i = 2; i < max_prime; i++)
	{
		if (is_prime[i])
		{
			primes[++primes[0]] = i;

			for (j = (i << 1); j < max_prime; j += i)
			{
				is_prime[j] = 0;
			}
		}
	}

	primes[primes[0] + 1] = 0x7FFFFFFF;
}


void regular_sieve(int max_prime, long blocksize, int ** factors, const int * primes, const int flags)
{
	long i, j;

	int * f;

	factors[0][0] = 1;
	factors[0][1] = 0;
	factors[1][0] = 1;
	factors[1][1] = 1;

	for (i = 0; i < blocksize; i++)
	{
		factors[i][0] = 0;
	}

	for (int i = 1, p = primes[i]; p < max_prime; p = primes[++i])
	{
		for (j = p; j < blocksize; j += p)
		{
			f = factors[j];
			f[++f[0]] = (flags & WITH_INDICES) ? i : p;
		}
	}
}


void segmented_sieve(int max_prime, long blocksize, long l, int ** factors, const int * primes, const int flags)
{
	if (l == 0)
	{
		regular_sieve(max_prime, blocksize, factors, primes, flags);
		return;
	}

	int j, k, offset, * f;

	for (k = 0; k < blocksize; k++)
	{
		factors[k][0] = 0;
	}

	for (int i = 1, p = primes[i]; p < max_prime; p = primes[++i])
	{
		offset = (p - (l % p)) % p;

		for (j = offset; j < blocksize; j += p)
		{
			f = factors[j];
			f[++f[0]] = (flags & WITH_INDICES) ? i : p;
		}
	}
}


void mod_sieve(const long blocksize, const long l, int ** factors, const int * primes, const int a, const int m)
{
	int i = 2, j, * f;

	long offset;
	long init_offset = ceil(((double) (l - a)) / m);

	//long total = floor(((double) (l + blocksize - a)) / m) - init_offset + 1;

	for (long k = 0; k < blocksize; k++)
	{
		factors[k][0] = 0;
	}

	const long max_prime = (long) sqrt(l + blocksize * m);
	long p = primes[i++];

	while (p < max_prime)
	{
		offset = (crt(0, p, a, m, l) - a) / m - init_offset;

		for (j = offset; j < blocksize; j += p)
		{
			f = factors[j];
			f[++f[0]] = p;
		}

		p = primes[i++];
	}
}

