/*
 * functions.h
 *
 *  Created on: Apr 30, 2014
 *      Author: antonmosunov
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#define WITH_INDICES 1

#ifndef ABS
#define ABS(X) (((X) < 0) ? (-(X)) : (X))
#endif

#ifndef MIN
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#endif


// HASH TABLE DECLARATIONS

typedef struct
{
	void * item;
	void * next;
} hentry_t;

typedef struct
{
	// hash_table.hpp members
	size_t item_size;					// size of item in hentry_t
	size_t size;						// number of buckets in the hash table
	size_t cur_size;					// current number of elemets
	hentry_t ** buckets;					// array of data buckets
	void * last_one;					// pointer to most recent entry

	// indexed_hash_table.hpp members
    long allocated;						// size of array
    hentry_t ** IDX;						// array of pointers to list elements
    long outstyle;						// output style (0=list, 1=hash table)

	int (*hash) (const void *);				// hash function to use
	int (*eq) (const void *, const void *);	// function for element comparison
	void (*del) (void *);					// called when item gets deleted
} htab_t;

void htab_init(htab_t * t, const size_t prime, const size_t item_size, int (*hash)(const void *), int (*eq) (const void *, const void *), void (*del) (void *));

void htab_empty(htab_t * t);

void htab_clear(htab_t * t);

void htab_insert(htab_t * t, const void * G);

void * htab_find(const htab_t * t, const void * G);

void htab_delete_from(htab_t * t, const long i);

void htab_delete(htab_t * t, const void * G);

void * htab_get(const htab_t * t, const long i);

size_t htab_size(const htab_t * t);


// VECTOR DECLARATIONS

typedef struct
{
	size_t alloc;
	size_t size;
	int * v;
} vec_t;

void vec_init(vec_t * v);

void vec_push_back(vec_t * v, int val);

void vec_set(vec_t * v, int val, const size_t j);

void vec_copy(vec_t * to, const vec_t * from);

void vec_copy_arr(vec_t * to, const int * from, const size_t size);

void vec_add(int * r, const vec_t * v, const vec_t * w);

void vec_zero(vec_t * v);

void vec_print(const vec_t * v);

void vec_clear(vec_t * v);


// MATRIX DECLARATIONS

typedef int ** mat_t;

mat_t mat_init(const size_t size);

void mat_print(mat_t M, const size_t size);

void mat_clear(mat_t M, const size_t size);

void smith_normal_form(mat_t M, const size_t size);


// MISCELLANEOUS FUNCTIONS

int next_prime(const int n);

// computes x = a (mod m), x = b (mod n), which satisfies x >= min
long crt(const int a, const int m, const int b, const int n, const long min);

char kronecker_symbol(long a, long p);

long divisors_list(long * result, long D, const int * pfactors, const int * primes, const int flags);

#endif /* FUNCTIONS_H_ */
