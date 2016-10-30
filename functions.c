/*
 * functions.c
 *
 *  Created on: Apr 30, 2014vec
 *      Author: antonmosunov
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <liboptarith/primes.h>
#include <liboptarith/gcd/gcd_binary_l2r.h>
#include "functions.h"

#ifndef ABS
#define ABS(X) ((X < 0) ? -X : X)
#endif




// HASH TABLE IMPLEMENTATION

void htab_init(htab_t * t, const size_t prime, const size_t item_size, int (*hash) (const void *), int (*eq) (const void *, const void *), void (*del) (void *))
{
	t->size = prime;
	t->item_size = item_size;
	t->cur_size = 0;
	t->buckets = (hentry_t **) calloc(prime, sizeof(void *));
	t->last_one = NULL;

	t->IDX = (hentry_t **) malloc(sizeof(hentry_t *));
	t->allocated = 1;
	t->outstyle = 0;

	t->hash = hash;
	t->eq = eq;
	t->del = del;
}

// hash_table<T>::empty()
void htab_empty(htab_t * t)
{
	register long i, end;

	end = t->cur_size;
	for (i = end - 1; i >= 0; i--)
	{
		htab_delete_from(t, i);
	}
}

// hash_table<T>::empty()
void htab_clear(htab_t * t)
{
	htab_empty(t);
	free(t->IDX);
	free(t->buckets);
}


// hash_table<T>::hash(const T & G)
void htab_insert(htab_t * t, const void * G)
{
	register long i;
	hentry_t * ptr, * nptr, * newone;

	// allocate new list element and a copy of G
	newone = (hentry_t *) malloc(sizeof(hentry_t));
	newone->item = (void *) malloc(t->item_size);
	memcpy(newone->item, G, t->item_size);
	newone->next = NULL;

	// compute correct bucket and append G to the end of the linked list
	i = t->hash(G) % t->size;

	if (i < 0)
	{
		i += t->size;
	}

	ptr = t->buckets[i];

	if (ptr)
	{
		while(ptr)
		{
			nptr = ptr;
			ptr = nptr->next;
		}
		nptr->next = newone;
	}
	else
	{
		t->buckets[i] = newone;
	}

	// append the pointer to the list of pointers
	if (t->cur_size == t->allocated)
	{
		// allocated more storage for the sequential list
		t->allocated <<= 1;
		t->IDX = (hentry_t **) realloc(t->IDX, t->allocated * sizeof(hentry_t *));
	}

	t->IDX[t->cur_size++] = newone;
}


void * htab_find(const htab_t * t, const void * G)
{
	void * target, * tptr;
	hentry_t * ptr;

	// compute correct bucket and perform linked list search
	long i = t->hash(G) % t->size;

	if (i < 0)
	{
		i += t->size;
	}

	target = NULL;

	ptr = t->buckets[i];

	while (ptr && !target)
	{
		tptr = ptr->item;

		if (t->eq(G, tptr))
		{
			target = tptr;
		}
		else
		{
			ptr = ptr->next;
		}
	}

	return target;
}


void htab_delete_from(htab_t * t, const long i)
{
	hentry_t * ptr, * pptr, * nptr, * temp;
	register long j;

	ptr = t->IDX[i];

	for (j = i; j < t->cur_size - 1; j++)
	{
		t->IDX[j] = t->IDX[j + 1];
	}
	t->IDX[t->cur_size - 1] = NULL;

	// find previous pointer in linked list
	j = t->hash(ptr->item) % t->size;

	if (j < 0)
	{
		j += t->size;
	}

	temp = t->buckets[j];
	pptr = NULL;

	while (temp != ptr)
	{
		pptr = temp;
		temp = temp->next;
	}

	// delete
	nptr = ptr->next;
	if (pptr)
	{
		pptr->next = nptr;
	}
	else
	{
		// G was the first element in the list - handle specially
		t->buckets[j] = nptr;
	}

	t->del(ptr->item);
	free(ptr->item);
	free(ptr);
	t->cur_size--;
}


void htab_delete(htab_t * t, const void * G)
{
	register long j;
	void * target;

	// check whether G is in the table
	target = htab_find(t, G);

	// if so, delete it
	if (target)
	{
		// compute j, its sequential index
		for (j = 0; j < t->cur_size; j++)
		{
			if (t->eq(t->IDX[j]->item, target))
			{
				break;
			}
		}

		// remove element j
		htab_delete_from(t, j);
	}
}

void * htab_get(const htab_t * t, const long i)
{
	return t->IDX[i]->item;
}

size_t htab_size(const htab_t * t)
{
	return t->cur_size;
}




// VECTOR IMPLEMENTATION

void vec_init(vec_t * v)
{
	v->alloc = 1;
	v->size = 0;
	v->v = (int *) malloc(sizeof(int));
}

void vec_clear(vec_t * v)
{
	free(v->v);
}

void vec_push_back(vec_t * v, int val)
{
	if (v->size == v->alloc)
	{
		v->alloc <<= 1;		// double the vector size
		v->v = (int *) realloc(v->v, v->alloc * sizeof(int));
	}

	v->v[v->size++] = val;
}

void vec_zero(vec_t * v)
{
	v->size = 0;
}

void vec_set(vec_t * v, int val, const size_t j)
{
	if (j + 1 > v->alloc)
	{
		v->alloc = j + 1;
		v->v = (int *) realloc(v->v, v->alloc * sizeof(int));
	}

	if (j + 1 > v->size)
	{
		memset(v->v + v->size, 0, (j + 1 - v->size) * sizeof(int));
		v->size = j + 1;
	}

	v->v[j] = val;
}

void vec_copy(vec_t * to, const vec_t * from)
{
	if (to->alloc < from->size)
	{
		to->alloc = from->alloc;
		to->v = (int *) realloc(to->v, to->alloc * sizeof(int));
	}

	to->size = from->size;
	memcpy(to->v, from->v, to->size * sizeof(int));
}

void vec_copy_arr(vec_t * to, const int * from, const size_t size)
{
	if (to->alloc < size)
	{
		to->alloc = size;
		to->v = (int *) realloc(to->v, to->alloc * sizeof(int));
	}

	to->size = size;
	memcpy(to->v, from, to->size * sizeof(int));
}

void vec_add(int * r, const vec_t * v, const vec_t * w)
{
	register size_t k = 0, ws = w->size, vs = v->size;
	if (ws < vs)
	{
		for (; k < ws; r[k] = v->v[k] + w->v[k], k++);
		for (; k < vs; r[k] = v->v[k], k++);
	}
	else
	{
		for (; k < vs; r[k] = v->v[k] + w->v[k], k++);
		for (; k < ws; r[k] = w->v[k], k++);
	}
}

void vec_print(const vec_t * v)
{
	printf("%d", v->v[0]);
	for (size_t i = 1; i < v->size; i++)
	{
		printf(" %d", v->v[i]);
	}
}




// MATRIX IMPLEMENTATION

mat_t mat_init(const size_t size)
{
	mat_t M = (mat_t) malloc(size * sizeof(int *));

	for (size_t i = 0; i < size; i++)
	{
		M[i] = (int *) calloc(size, sizeof(int));
	}

	return M;
}


void mat_print(mat_t M, const size_t size)
{
	for (size_t i = 0; i < size; i++)
	{
		printf("%d", M[i][0]);
		for (long j = 1; j < size; j++)
		{
			printf("\t%d", M[i][j]);
		}
		printf("\n");
	}
}


void mat_clear(mat_t M, const size_t size)
{
	for (size_t i = 0; i < size; i++)
	{
		free(M[i]);
	}

	free(M);
}


void smith_normal_form(mat_t M, const size_t size)
{
	int i, c;
	int j, k, l;			// looping variables
	int b;					// a varible used to test when we are done
	int64_t dd = 0, u, v;	// values used when computing the extended gcd
	int r;					// used as a remainder after division
	char done = 0;			// to signify completion or not

	if (size == 1)
	{
		return;
	}

	long n = size;
	long temp;			// a temporary vector to hold columns and rows
	i = n - 1;

	//while (done == 0)
	for (i = n - 1; i > 0; i = (done == 1) ? i - 1 : i)
	{
		do
		{
			//initialize j and c for row reduction
			c = 0;
			for (j = i - 1; j >= 0; j--)
			{
				if (M[i][j] == 0)
				{
					continue;
				}

				//start step 4 of the algorithm
				if (M[i][j] % M[i][i] == 0)
				{
					u = 1;
					v = 0;
					dd = M[i][i];
				}
				else
				{
					dd = xgcd_binary_l2r_s64(&u, &v, M[i][i], M[i][j]);
				}

				// use r and b as temporary variables to avoid
				// excessive division in the for loop.
				r = M[i][i] / dd;
				b = M[i][j] / dd;

				for (k = 0; k < n; k++)
				{
					temp = u * M[k][i] + v * M[k][j];
					M[k][j] = r * M[k][j] - b * M[k][i];
					M[k][i] = temp;
				}
			}
#ifdef DEBUG
printf("Done Transformation\n");
printf("dd = %ld, u = %ld, v = %ld\nM=\n", dd, u, v);
mat_print(M, size);
#endif
			//initialize j and c for column reduction
			for (j = i - 1; j >= 0; j--)
			{
				if (M[j][i] == 0)
				{
					continue;
				}

				//start step 7 of the algorithm
				if (M[j][i] % M[i][i] == 0)
				{
					u = 1;
					v = 0;
					dd = M[i][i];
				}
				else
				{
					dd = xgcd_binary_l2r_s64(&u, &v, M[i][i], M[j][i]);
				}

				r = M[i][i] / dd;
				b = M[j][i] / dd;

				for (k = 0; k < n; k++)
				{
					temp = u * M[i][k] + v * M[j][k];
					M[j][k] = r * M[j][k] - b * M[i][k];
					M[i][k] = temp;
				}

#ifdef DEBUG
printf("Done Transform 2\n");
printf("dd = %ld, u = %ld, v = %ld\nM=\n", dd, u, v);
mat_print(M, size);
#endif
				c++;
			}
		} while (c > 0);

		//now do step 9
		b = M[i][i];
		done = 1;

		for (k = 0; k < i && done == 1; k++)
		{
			for (l = 0; l < i && done == 1; l++)
			{
				r = M[k][l] % b;
				if (r != 0)
				{
					for (long t = 0; t < n; t++)
					{
						M[i][t] += M[k][t];
					}
					done = 0;
				}
			}
		}

#ifdef DEBUG
printf("Done Step 9\nM=\n");
mat_print(M, size);
#endif
	}
}


int next_prime(const int n)
{
	for (int i = 0; i < prime_list_count; i++)
	{
		if (prime_list[i] > n)
		{
			return prime_list[i];
		}
	}

	return -1;
}


// computes x = a (mod m), x = b (mod n), which satisfies x >= min
// considers only the case a = 0, n = 2^k
long crt(const int a, const int m, const int b, const int n, const long min)
{
	int c, d;

	xgcd_binary_l2r_s32(&c, &d, m, n);

	long x = a + (b - a) * c * m;

	long k = (long) ceil(((double) (min - x)) / (m * n));

	return (x + m * n * k);
}


char kronecker_symbol(long a, long p)
{
	if (p == 2)
	{
		if (a & 1)
		{
			return ((a & 7) == 1 || (a & 7) == 7) ? 1 : -1;
		}
		else
		{
			return 0;
		}
	}

	char is_negative = 0;
	if (a < 0)
	{
		is_negative = ((p & 3) == 3);
		a = -a;
	}

	long temp;
	char t = 1;

	while (a != 0)
	{
		while ((a & 1) == 0)
		{
			a >>= 1;
			if ((p & 7) == 3 || (p & 7) == 5)
			{
				t = -t;
			}
		}

		if (a < p)
		{
			temp = a;
			a = p;
			p = temp;

			if ((a & 3) == 3 && (p & 3) == 3)
			{
				t = -t;
			}
		}

		a = ((a - p) >> 1);

		if ((p & 7) == 3 || (p & 7) == 5)
		{
			t = -t;
		}
	}

	if (is_negative)
	{
		t = -t;
	}

	return (p == 1) ? t : 0;
}


long divisors_list(long * result, long D, const int * pfactors, const int * primes, const int flags)
{
	if (D == 0)
	{
		return 0;
	}

	if (D == 1)
	{
		result[0] = 1;
		return 1;
	}

	long total = pfactors[0];
	long i, j, p;

	long prev = 0, cur = 1;
	long p_pow;
	result[0] = 1;

	for (i = 1; i <= total; i++)
	{
		prev = cur;
		p = (flags & WITH_INDICES) ? primes[pfactors[i]] : pfactors[i];
		p_pow = 1;
		while (D % p == 0)
		{
			D /= p;
			p_pow *= p;
			for (j = 0; j < prev; j++)
			{
				result[cur++] = result[j] * p_pow;
			}
		}
	}

	if (D > 1)
	{
		prev = cur;
		for (j = 0; j < prev; j++)
		{
			result[cur++] = result[j] * D;
		}
	}

	return cur;
}
