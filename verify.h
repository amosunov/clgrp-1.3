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

#ifndef VERIFY_H_
#define VERIFY_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <gmp.h>
//#include <omp.h>

#ifdef WITH_PARI
#include <pari/pari.h>
#endif

#include "functions.h"
#include "sieve.h"

// corresponds to 963761198400, smallest # < 2^40 which has largest # of divisors
#define MAX_DIVISORS 12000 //6720

void partial_left_hand_side(mpz_t LHS, const char * file, const char * folder, const int index, const long blocksize, const long x, const int a, const int m, const int * primes, int ** factors);

void left_hand_side(mpz_t LHS, const long D_max, const long files, const int * primes, const char * folder);

void partial_right_hand_side(mpz_t sum, const long blocksize, const long l, const int * primes, int ** factors);

void right_hand_side(mpz_t RHS, const long n_max, const long blocksize, const int * primes);

#endif /* SIEVE_H_ */
