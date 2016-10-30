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

#ifndef CLASS_GROUP_H_
#define CLASS_GROUP_H_

#include <libqform/s64_qform.h>

#include "functions.h"

#define MAX_RANK 10


// FORM

typedef struct
{
	s64_qform_t form;
	vec_t value;		// for hash table
} form_t;

static __inline__
int hash_form_t(const void * a)
{
	return ((form_t *) a)->form.a;
}

static __inline__
int eq_form_t(const void * a, const void * b)
{
	s64_qform_t * t = (s64_qform_t *) a;
	s64_qform_t * s = &((form_t *) b)->form;
	return ((t->a == s->a) && (t->b == s->b));
}

static __inline__
void del_form_t(void * a)
{
	vec_clear(&((form_t *) a)->value);
}


// CLASS GROUP COMPUTATION / TABULATION

#ifdef WITH_PARI
#include <pari/pari.h>

void pari_verify(int * result, const long D);
#endif

int next(group_pow_t * gp, form_t * R, const int init_pow, int prime_index);

int h_upper_bound(const long D);

int h_lower_bound(const long D);

int compute_group_bjt(int * result, const long D, const int init_pow, const int h_star, htab_t * R, htab_t * Q);

void tabulate_bjt(const int index, const long D_total, const char * file, const char * folder, const int a, const int m,
				const int * small_primes, int ** h_factors, int ** D_factors, int * h_list);

#endif /* CLASS_GROUP_H_ */
