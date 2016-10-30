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
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include <liboptarith/primes.h>
#include <liboptarith/sqrtmodp_list.h>

#include "clgrp.h"
#include "sieve.h"

#ifdef WITH_PARI
#include <pari/pari.h>
#endif

#ifndef MIN
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#endif

#ifndef MAX
#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#endif

#ifndef ABS
#define ABS(X) (((X) < 0) ? (-(X)) : (X))
#endif

#define ADDR (form_t *) malloc(sizeof(form_t))

// from /ANTL/include/ANTL/quadratic/quadratic_order.hpp, table OQvals[]
const long Q_table[] =
{	947, 2269, 3929, 6011, 8447, 11093, 14149, 17393, 20921, 24733, 28807,
	33151, 37619, 42533, 47507, 52859, 58321, 64231, 70099, 76463 };
  
const double C_table[] =
{
	9785883.57965490035712718963623046875,
	62938341.824559591710567474365234375,
	201442194.1823149621486663818359375,
	494557498.147336781024932861328125,
	1013054914.3328609466552734375
};


#ifdef WITH_PARI
void pari_verify(int * result, const long D)
{
	GEN g = stoi(D);
	GEN v = quadclassunit0(g, 0, 0, 0);
	GEN clgrp = compo(v, 2);

	for (int i = 1; i < lg(clgrp); i++)
	{
		if (result[i - 1] != gtolong(compo(clgrp, i)))
		{
			char err[100];
			sprintf(err, "D=%ld, Cl[%d] = %d, actual %ld\n", D, i - 1, result[i - 1], gtolong(compo(clgrp, i)));
			perror(err);
			exit(1);
		}
	}

	cgiv(clgrp);
	cgiv(v);
	cgiv(g);
}
#endif

void qform_pow_s32(group_pow_t * pow, qform_t * R, const qform_t * A, int32_t exp)
{
	if (exp < 0)
	{
		s64_qform_group_t * group = (s64_qform_group_t *) pow->group;
		s64_qform_t B;
		s64_qform_set(group, &B, A);
		s64_qform_inverse(group, &B);
		qform_pow_u32(pow, R, &B, -exp);
	}
	else
	{
		qform_pow_u32(pow, R, A, exp);
	}
}

	
// init_pow = class_number / order of the p-group
int next(group_pow_t * gp, form_t * R, const int init_pow, int prime_index)
{
	s64_qform_group_t * group = (s64_qform_group_t *) gp->group;

	s64_qform_set_id(group, &R->form);

	while (s64_qform_is_id(group, &R->form))
	{
		while (s64_qform_is_primeform(group, &R->form, prime_list[++prime_index]) == 0);

		qform_pow_s32(gp, &R->form, &R->form, init_pow);
	}

	return prime_index;
}

// uses Olivier Ramare's upper bounds on L(1,x)
int h_upper_bound(const long D)
{
	double E;
	if ((D & 3) == 0)
	{
		E = 0.25 * log(ABS(D)) + 1.25 - 0.5 * log(3);
	}
	else
	{
		E = 0.5 * log(ABS(D)) + 2.5 - log(6);
	}
	
	return (exp(E) * sqrt(ABS(D))) / M_PI;
}

// Algorithm from ANTL/src/L_function/L_function_long.cpp,
// approximateL1_impl(2) routine
int h_lower_bound(const long D)
{
	long Q, Q2, P, kron, r, prime_index;
	double C, E, wt;
	register long i;
	
	const long D_abs = ABS(D);

	r = ((long) (log10(D_abs)) / 5);
	C = C_table[r];
    Q = Q_table[r];
	Q2 = Q << 1;			// Q2 = Q*2;

    // compute partial product  p < Q
    E = 0.0;
    prime_index = 0;
    P = prime_list[prime_index];		// primes start at 2
    while (P < Q)
    {
        kron = ((r = D % P) != 0) ? ((sqrtmodp[P][r + P] > 0) ? 1 : -1) : 0;
        E += log(((double) P) / (P - kron));
        P = prime_list[++prime_index];
    }

    // computed weighted partial products for Q < p < 2Q
    wt = 1.0;
    for (i = Q; i <= P; i++)
    {
		wt -= i * log(i) / C;
	}

    while (P < Q2)
	{
		kron = ((r = D % P) != 0) ? ((sqrtmodp[P][r + P] > 0) ? 1 : -1) : 0;
		E += wt * log(((double) P) / (P - kron));
		P = prime_list[++prime_index];
		wt -= ((P - 1) * log(P - 1) + P * log(P)) / C;
	}
	
	double h = (exp(E) * sqrt(D_abs)) / (M_SQRT2 * M_PI);
	
	if ((D_abs & 7) == 3)
	{
		h /= 3.0;
	}
	
	int h_star = (int) h;

	if (h_star < 5)
	{
		h_star++;
	}

	return h_star;
}


// Ramachandran's thesis, p. 47
int compute_group_bjt(int * result, const long D, const int init_pow, const int h_star, htab_t * R, htab_t * Q)
{
#ifdef DEBUG
printf("h_star=%d\n", h_star);
#endif

	mat_t M = mat_init(MAX_RANK);

	s64_qform_group_t group;
	s64_qform_group_init(&group);
	s64_qform_group_set_discriminant_s64(&group, D);

	group_pow_t gp;
	group_pow_init(&gp, &group.desc.group);

	form_t * e, * it;

	form_t ne;
	s64_qform_set_id(&group, &ne.form);
	vec_init(&ne.value);
	vec_push_back(&ne.value, 0);
	htab_insert(R, &ne);
	vec_init(&ne.value);
	vec_push_back(&ne.value, 0);
	htab_insert(Q, &ne);

	int omega = 2, h = 1, det = 1, s, y, u, i, j = 0, k, prime_index = -1, rank = 0, q, n, t;
	size_t R_size, R_prev_size, Q_size, cur_index;

	s64_qform_t g, g_inv, a, b, c, temp;

	char is_break;

	while (h < h_star)
	{
#ifdef DEBUG
printf("\n\n\nNEW ITERATION: j=%d\n", j);
printf("INITIAL b:\n");
for (int p = 0; p <= j; p++)
{
	printf("b_vec[%d]=", p);
	for (int l = 0; l <= p; l++)
	{
		printf(" %d", M[p][l]);
	}
	printf("\n");
}
#endif

		Q_size = htab_size(Q);
		R_prev_size = htab_size(R);
		cur_index = R_prev_size;

		// initializations
		prime_index = next(&gp, &ne, init_pow, prime_index);
		s64_qform_set(&group, &g, &ne.form);

#ifdef DEBUG
printf("g[%d]=", j);
s64_qform_print(&group, &g);
printf(", <- prime_index=%d\n", prime_index);
#endif

		s = 1;
		y = omega;
		u = omega;
		s64_qform_set(&group, &g_inv, &g);
		s64_qform_inverse(&group, &g_inv);
		s64_qform_set_id(&group, &a);
		qform_pow_u32(&gp, &b, &g, omega);
		s64_qform_set(&group, &c, &b);

		is_break = 0;

		// check if current generator is contained in current subgroup
		if (j > 0)
		{
#ifdef DEBUG
printf("\n\nContained in current subgroup?\n");
#endif

			for (i = 0; i < Q_size; i++)
			{
				is_break = 0;
				it = (form_t *) htab_get(Q, i);
				
				s64_qform_compose(&group, &temp, &it->form, &g);
				e = (form_t *) htab_find(R, &temp);

				if (e != NULL)
				{
					n = MIN(MIN(e->value.size, it->value.size), j);
					for (k = 0; k < n; k++)
					{
						if (it->value.v[k] + e->value.v[k] >= M[k][k])
						{
							is_break = 1;
							break;
						}
					}

					if (!is_break)
					{
#ifdef DEBUG
printf("Yes! ");
s64_qform_print(&group, &e->form);
printf("\n");
#endif

						M[j][j] = 1;
						break;
					}
				}
			}
		}

		while (M[j][j] == 0)
		{
#ifdef DEBUG
printf("\nBaby steps, s=%d, u=%d\n", s, u);
#endif
			// compute new baby steps
			for (i = s; i <= u; i++)
			{
				s64_qform_compose(&group, &a, &a, &g_inv);

				if (a.a == 1)
				{
					M[j][j] = i;
					break;
				}

				if (s == 1 && i > 1)
				{
					for (t = 0; t < R_prev_size; t++)
					{
						it = (form_t *) htab_get(R, t);
				
#ifdef DEBUG
printf("From hash table we got form ");
s64_qform_print(&group, &it->form);
printf("\n");
#endif
						is_break = 0;
						s64_qform_compose(&group, &temp, &it->form, &a);
						e = (form_t *) htab_find(Q, &temp);

						if (e != NULL)
						{
							n = MIN(MIN(e->value.size, it->value.size), j);
							for (k = 0; k < n; k++)
							{
								if (it->value.v[k] + e->value.v[k] >= M[k][k])
								{
									is_break = 1;
									break;
								}
							}

							if (!is_break)
							{
#ifdef DEBUG
printf("s=%d, i=%d, Found ", s, i);
s64_qform_print(&group, &e->form);
printf(" = ");
s64_qform_print(&group, &it->form);
printf(" * ");
s64_qform_print(&group, &a);
printf("\n");
#endif
								vec_add(M[j], &it->value, &e->value); 
								M[j][j] = i;
								
#ifdef DEBUG
printf("b_vec[%d]=", j);
for (int l = 0; l <= j; l++)
{
printf(" %d", M[j][l]);
}
printf("\n");
#endif
								goto giant_steps;
							}
						}
						else
						{
							is_break = 1;
						}

						if (is_break)
						{
							s64_qform_set(&group,&ne.form, &temp);
							vec_init(&ne.value);
							vec_copy(&ne.value, &it->value);
							vec_set(&ne.value, i, j);
							htab_insert(R, &ne);
							cur_index++;

#ifdef DEBUG
printf("s=%d, i=%d, Adding ", s, i);
s64_qform_print(&group, &temp);
printf(" <->");
for (int l = 0; l <= j; l++)
{
printf(" %d", ne.value.v[l]);
}
printf(" to R\n");
#endif
						}
					}
				}
				else
				{
					for (t = 0; t < R_prev_size; t++)
					{
						it = (form_t *) htab_get(R, t);
				
#ifdef DEBUG
printf("t=%d, R_prev_size=%lu, From hash table we got form ", t, R_prev_size);
s64_qform_print(&group, &it->form);
printf("\n");
#endif
						s64_qform_compose(&group, &ne.form, &it->form, &a);
						vec_init(&ne.value);
						vec_copy(&ne.value, &it->value);
						vec_set(&ne.value, i, j);

#ifdef DEBUG
printf("s=%d, i=%d, Adding ", s, i);
s64_qform_print(&group, &ne.form);
printf(" <->");
for (int l = 0; l <= j; l++)
{
printf(" %d", ne.value.v[l]);
}
printf(" to R\n");
#endif

						htab_insert(R, &ne);
						cur_index++;
					}
				}
			}

			// compute giant steps
			giant_steps:
#ifdef DEBUG
printf("\nGiant steps\n");
#endif
			while (M[j][j] == 0 && y < u * u)
			{
#ifdef DEBUG
printf("y=%d, u=%d, b = ", y, u);
s64_qform_print(&group, &b);
printf("\n");
#endif

				for (t = 0; t < Q_size; t++)
				{
					it = htab_get(Q, t);
					is_break = 0;
					s64_qform_compose(&group, &temp, &it->form, &b);

#ifdef DEBUG
printf("Looking for ");
s64_qform_print(&group, &temp);
printf(" = ");
s64_qform_print(&group, &it->form);
printf(" * ");
s64_qform_print(&group, &b);
printf("\n");
#endif					

					e = htab_find(R, &temp);

					if (e != NULL)
					{
						n = MIN(MIN(it->value.size, e->value.size), j);
						for (k = 0; k < n; k++)
						{
							if (it->value.v[k] + e->value.v[k] >= M[k][k])
							{
								is_break = 1;
								break;
							}
						}

						if (!is_break)
						{
							vec_add(M[j], &it->value, &e->value);
							M[j][j] += y;

							if (M[j][j] == 0)
							{
								printf("b_vec is zero!\n");
								continue;
							}

#ifdef DEBUG
printf("Found ");
s64_qform_print(&group, &e->form);
printf(" = ");
s64_qform_print(&group, &it->form);
printf(" * ");
s64_qform_print(&group, &b);
printf("\nb_vec[%d]=", j);
for (int l = 0; l <= j; l++)
{
	printf(" %d", M[j][l]);
}
printf("\n");
#endif
							goto double_step;
						}
					}
				}

				y += u;
				s64_qform_compose(&group, &b, &b, &c);
			}

			// double step width
			double_step:
#ifdef DEBUG
printf("Doubling step width...\n");
#endif
			s = u + 1;
			u *= 2;
			s64_qform_square(&group, &c, &c);
		}

		if (M[j][j] > 1)
		{

#ifdef DEBUG
printf("RESULT:\n");
for (int p = 0; p <= j; p++)
{
	printf("b_vec[%d]=", p);
	for (int l = 0; l <= p; l++)
	{
		printf(" %d", M[p][l]);
	}
	printf("\n");
}
#endif

			h *= M[j][j];

			if (h < h_star)
			{
				// update R and Q
				q = (int) ceil(sqrt(M[j][j]));
				det *= q;
				R_size = htab_size(R);

#ifdef DEBUG
printf("\n\nUpdating R and Q\nq=%d, det=%d, R_size=%lu\n", q, det, R_size);
#endif

				for (t = R_size - 1; t >= det; t--)
				{
#ifdef DEBUG
form_t * foo = (form_t *) htab_get(R, t);
printf("Removing ");
s64_qform_print(&group, &foo->form);
printf(" <-> ");
vec_print(&foo->value);
printf("\n");
#endif
					htab_delete_from(R, t);
				}

				Q_size = htab_size(Q);
				qform_pow_u32(&gp, &c, &g, q);
				s64_qform_set_id(&group, &b);

				for (i = 1; (i < q) && (i * q < M[j][j]); i++)
				{
					s64_qform_compose(&group, &b, &b, &c);
#ifdef DEBUG
printf("g^%d = ", i * q);
s64_qform_print(&group, &b);
printf("\n");
#endif
					for (t = 0; t < Q_size; t++)
					{
						it = (form_t *) htab_get(Q, t);
						s64_qform_compose(&group, &ne.form, &b, &it->form);
						vec_init(&ne.value);
						vec_copy(&ne.value, &it->value);
						vec_set(&ne.value, i * q, j);
						htab_insert(Q, &ne);

#ifdef DEBUG
printf("Adding ");
s64_qform_print(&group, &ne.form);
printf(" = ");
s64_qform_print(&group, &it->form);
printf(" * g^%d <->", i * q);
for (int l = 0; l <= j; l++)
{
	printf(" %d", ne.value.v[l]);
}
printf(" to Q\n");
#endif
					}			
				}
			}

			j++;
		}
		else
		{
			M[j][j] = 0;
		}
	}

	if (h > 2 * h_star)
	{
		printf("BUG! h > 2 * h_star, TERMINATING: D=%ld\n", D);
		fflush(stdout);
		exit(1);
	}

	if (h != 1)
	{
		rank = j;

		group_pow_clear(&gp);
		s64_qform_group_clear(&group);

#ifdef DEBUG
printf("ANSWER:\n");
for (int i = 0; i < rank; i++)
{
	for (int j = 0; j <= i; j++)
	{
		printf("%d ", M[i][j]);
	}
	printf("\n");
}
printf("Matrix before:\n");
mat_print(M, rank);
#endif

		smith_normal_form(M, rank);

#ifdef DEBUG
printf("Matrix after:\n");
mat_print(M, rank);
#endif

		j = 1;

		for (int i = 0; i < rank; i++)
		{
			if (M[i][i] > 1)
			{
				result[j++] = M[i][i];
			}
		}
	}
	else
	{
		result[1] = 1;
		j = 2;
	}

	htab_empty(R);
	htab_empty(Q);
	
	result[0] = h;

	mat_clear(M, MAX_RANK);

	return (j - 1);
}




void tabulate_bjt(const int index, const long D_total,
		const char * file, const char * folder,
		const int a, const int m,
		const int * primes, int ** h_factors,
		int ** D_factors, int * h_list)
{
	char name[500], data[200];
	int fd;

	struct timeval begin, end;
	unsigned long exec_time;

	sprintf(name, "%s/cl%dmod%d/cl%dmod%d.%d.gz", folder, a, m, a, m, index);
	if (access(name, F_OK) != -1)
	{
		printf("The file %s exists, thus terminating.\n", name);
		return;
	}

	if (file)
	{
		sprintf(name, "%s/%s%d", folder, file, index);

		if (index == 0 && (a == 4 || a == 3))
		{
			fd = open(name, O_RDWR);

			if (fd == -1)
			{
				sprintf(data, "Unable to open the file %s\n", name);
				perror(data);
				fflush(stderr);
				exit(1);
			}

			int * h_temp = (int *) mmap(0, sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

			if (h_temp == MAP_FAILED)
			{
				sprintf(data, "Unable to map the file %s for writing\n", name);
				perror(data);
				fflush(stderr);
				exit(1);
			}

			if (a == 4)
			{
				h_temp[0] = 2;
			}
			else
			{
				h_temp[1] = 3;
			}

			munmap(h_temp, sizeof(int));
			close(fd);
		}

		fd = open(name, O_RDONLY);

		if (fd == -1)
		{
			sprintf(data, "Unable to open the file %s for reading\n", name);
			perror(data);
			fflush(stderr);
			exit(1);
		}
	}
	else
	{
		fd = 0;
	}

	const long D_max = (index + 1) * D_total * m;
	long D_cur = 0, f_cur = FAC_TOTAL;
	int * D_cur_factors;
	long D_temp;

	int result[15];

	char is_discriminant = 1;

	int h, dist = 0, r, rank;

	// compute an upper bound on the size of the table
	int h_max = h_upper_bound(-D_max);

	int table_size = next_prime((((int) sqrt(h_max)) << 1) - 1);

	if (table_size == -1)
	{
		perror("Not enough primes in liboptarith/primes.h\n");
		fflush(stderr);
		exit(1);
	}

	htab_t R, Q;
	htab_init(&R, table_size, sizeof(form_t), &hash_form_t, &eq_form_t, &del_form_t);
	htab_init(&Q, table_size, sizeof(form_t), &hash_form_t, &eq_form_t, &del_form_t);
	sprintf(name, "%s/cl%dmod%d", folder, a, m);
	mkdir(name, 0744);
	sprintf(name, "%s/cl%dmod%d/cl%dmod%d.%d", folder, a, m, a, m, index);
	creat(name, 0744);
	FILE * clfd = fopen(name, "w");

	if (clfd == NULL)
	{
		sprintf(data, "Unable to open %s\n", name);
		perror(data);
		fflush(stderr);
		exit(1);
	}

	const long res = (((a & 3) != 3) ? a : 1);

	int init_pow = 1, h_fac, h_temp, h_fac_total;
	const int * h_cur_factors;

	gettimeofday(&begin, NULL);
	for (long D = index * D_total * m + a; D < D_max; D += m, D_cur++, f_cur++)
	{
		if (f_cur == FAC_TOTAL)
		{
			f_cur = 0;
			mod_sieve(FAC_TOTAL, D - a, D_factors, primes, a, m);
			read(fd, h_list, FAC_TOTAL * sizeof(int));
		}

		is_discriminant = 1;
		D_cur_factors = D_factors[f_cur];
		D_temp = D / res;

		for (int j = 1; j <= D_cur_factors[0]; j++)
		{
			D_temp /= D_cur_factors[j];
			if (D_temp % D_cur_factors[j] == 0)
			{
				is_discriminant = 0;
				break;
			}
		}

		if (is_discriminant)
		{
			if (file)
			{
				h = h_list[f_cur];

				if (a == 4)
				{
					h /= 2;
				}
				else if (a == 3)
				{
					h /= 3;
				}

				init_pow = 1;
				h_cur_factors = h_factors[h];
				h_fac_total = h_cur_factors[0];
				h_temp = h;

				for (int i = 1; i <= h_fac_total; i++)
				{
					h_fac = h_cur_factors[i];
					h_temp /= h_fac;
					if (h_temp % h_fac != 0)
					{
						init_pow *= h_fac;
					}
				}

				h /= init_pow;
			}
			else
			{
				h = h_lower_bound(-D);
			}

			//printf("%ld: D=%ld, #factors=%u, h=%d, init_pow=%d\n", D_cur, D, D_cur_factors[0], h, init_pow);
			//fflush(stdout);

			rank = compute_group_bjt(result, -D, init_pow, h, &R, &Q);

			h = result[0] * init_pow;
			result[1] *= init_pow;

			sprintf(name, "%d\t%u\t", dist, h);

			for (r = 1; r < rank; r++)
			{
				sprintf(data, "%d ", result[r]);
				strcat(name, data);
			}
			sprintf(data, "%d\n", result[rank]);
			strcat(name, data);

			#ifdef WITH_PARI
			pari_verify(result + 1, -D);
			#endif

			fputs(name, clfd);

			dist = 1;
		}
		else
		{
			dist++;
		}
	}

	fclose(clfd);
	sprintf(name, "gzip %s/cl%dmod%d/cl%dmod%d.%d", folder, a, m, a, m, index);
	system(name);

	htab_clear(&R);
	htab_clear(&Q);
	gettimeofday(&end, NULL);
	exec_time = (end.tv_sec * 1e6 + end.tv_usec) - (begin.tv_sec * 1e6 + begin.tv_usec);
	printf("index=%d, took %.3f\n", index, exec_time / 1e6);
	fflush(stdout);

	if (file)
	{
		close(fd);

		#ifndef KEEP_FILES
		sprintf(name, "%s/%s%d", folder, file, index);
		remove(name);
		#endif
	}
}
