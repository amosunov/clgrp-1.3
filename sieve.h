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

#ifndef SIEVE_H_
#define SIEVE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// 2^20=1048576
#define FAC_TOTAL 1048576

void prime_sieve(const int max_prime, int * primes);

void regular_sieve(const int max_prime, const long blocksize, int ** factors, const int * primes, const int flags);

void segmented_sieve(const int max_prime, const long blocksize, const long l, int ** factors, const int * primes, const int flags);

void mod_sieve(const long blocksize, const long l, int ** factors, const int * primes, const int a, const int m);

#endif /* SIEVE_H_ */
