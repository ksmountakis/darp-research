#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

typedef struct {
	int i;
	float e, l;
	int source;
	int sink;
} tuple_t;

typedef struct {
	int i;
	int leak;
	float s;
} chain_node_t;

typedef struct {
	int (*prioritize)(const void*, const void*);
	int *o;
	float *s;
	int *v;
	int *inserted;
	tuple_t *tuples;
	chain_node_t **chain;
} aux_t;

typedef struct {
	int N, n, m;
	float T, Q, L;
	int i, j, k, z;
	float *x, *y, *d, *q, *e, *l;
	double **t;
} darp_t;


int priority_criterion(const void *t1ptr, const void *t2ptr) {
	tuple_t *t1 = (tuple_t*)t1ptr;
	tuple_t *t2 = (tuple_t*)t2ptr;
	float time1 = t1->l;
	float time2 = t2->l;

	if (t1->source || t2->sink)
		return -1;

	if (t2->sink || t1->source)
		return +1;

	if (time1 < time2)
		return -1;

	if (time1 == time2)
		return 0;
	return +1;
}

static void
initialize(aux_t *aux, darp_t darp)
{
	int k;

	aux->prioritize = priority_criterion;
	assert((aux->tuples = calloc(sizeof(tuple_t), darp.N)));
	assert((aux->v = calloc(sizeof(int), darp.n)));
	assert((aux->inserted = calloc(sizeof(int), darp.N)));
	assert((aux->o = calloc(sizeof(int), darp.N)));
	assert((aux->chain = calloc(sizeof(chain_node_t*), darp.m)));
	for (k=0; k < darp.m; k++) 
		assert((aux->chain[k] = calloc(sizeof(chain_node_t), darp.N)));
}

static void
prioritize(aux_t *aux, const darp_t darp, float *tardvec, int *a)
{
	int i, z;
	// sort according to priority criterion
	for (i=0; i < darp.N; i++) {
		aux->tuples[i].i = i;
		aux->tuples[i].e = darp.e[i];
		aux->tuples[i].l = darp.l[i] - 0.5*tardvec[i]*rand()/RAND_MAX;
		aux->tuples[i].source = (i == 0);
		aux->tuples[i].sink = (i == darp.N-1);
	}
	qsort(aux->tuples, darp.N, sizeof(tuple_t), aux->prioritize);

	// initialize priority rule
	for (z=0; z < darp.N; z++) {
		a[z] = aux->tuples[z].i;
		aux->o[aux->tuples[z].i] = z;
	}
}

static void
rho(aux_t *aux, const darp_t darp, int *a, float *s, float ***arcmat, float *tardvec, float *tardmax, float *tardsum)
{
	int i, j, k, p, kmin, kmax;
	int num_inserted;
	float tard, best_tard;
	int best_k, best_found;

	// initialize chain-form schedule
	for (k=0; k < darp.m; k++) {
		aux->chain[k][0].i = 0;
		aux->chain[k][0].leak = darp.Q;
	}

	// generate schedule
	s[0] = 0;
	num_inserted = 0;
	memset(aux->inserted, 0, sizeof(int)*darp.N);
	memset(tardvec, 0, sizeof(float)*darp.N);

	*tardsum = 0;
	*tardmax = 0;
	while (num_inserted < darp.N-2)
	{
		// pick next eligible event from priority rule a[]
		for (p=0; p < darp.N; p++) {
			j = a[p];

			if (j == 0 || j == darp.N-1) 
				continue;

			//printf("\n%d ", j);

			if (aux->inserted[j]) {
				//printf("already inserted");
				continue;
			}

			if (j > darp.n && j < 2*darp.n+1) {
				if (!aux->inserted[j-darp.n]) {
					//printf("pickup not yet inserted");
					continue;
				}
				kmin = kmax = aux->v[j-darp.n];
			} else {
				kmin = 0;
				kmax = darp.m-1;
			}

			best_tard = 10*darp.T;
			best_found = best_k = 0;

			for (k=kmin; k <= kmax; k++) {
				if (aux->chain[k][0].leak < darp.q[j]) 
					continue;

				i = aux->chain[k][0].i;
				s[j] = s[i] + darp.t[i][j] + darp.d[i];
				s[j] = fmax(darp.e[j], s[j]);
				tard = s[j] - darp.l[j];

				if (tard <= best_tard) {
					best_k = k;
					best_tard = tard;
					best_found = 1;
				}
			}

			if (best_found) {
				*tardsum += darp.t[aux->chain[best_k][0].i][j];
				*tardmax = fmax(*tardmax, best_tard);
				//printf("(%d,%d):%d\n", chain[best_k][0].i, j, k);
				arcmat[aux->chain[best_k][0].i][j][best_k] = 1;
				aux->chain[best_k][0].i = j;
				aux->chain[best_k][0].leak -= darp.q[j];
				aux->inserted[j] = 1;
				tardvec[j] = best_tard;
				if (j <= darp.n)
					aux->v[j] = best_k;
				num_inserted++;
				//printf("inserted %d (k=%d leak=%d tard=%.1f)\n", j, best_k, chain[best_k][0].leak, best_tard);
			} else {
				//printf("no sufficient leak");
			}
		}
	}

	// append sink to all machines
	for (k=0; k < darp.m; k++) {
		i = aux->chain[k][0].i;
		arcmat[i][darp.N-1][k] = 1;
		s[darp.N-1] = fmax(s[darp.N-1], s[i] + darp.d[i] + darp.t[i][darp.N-1]);
	}
}

static void
darp_create(darp_t *darp)
{
	int i, j, err;
	// <m> <n> <max veh ride minutes Tk> <veh capacity Qk> <max user ride minutes Li>
	err = scanf("%d %d %f %f %f\n", &darp->m, &darp->N, &darp->T, &darp->Q, &darp->L);
	if (err != 5) {
		fprintf(stderr, "error processing line 1\n");
		exit(EXIT_FAILURE);
	}

	darp->N = 2*darp->N + 2;
	darp->n = (darp->N-2)/2;

	assert((darp->x = calloc(sizeof(float), darp->N)));
	assert((darp->y = calloc(sizeof(float), darp->N)));
	assert((darp->d = calloc(sizeof(float), darp->N)));
	assert((darp->q = calloc(sizeof(float), darp->N)));
	assert((darp->e = calloc(sizeof(float), darp->N)));
	assert((darp->l = calloc(sizeof(float), darp->N)));
	assert((darp->t = calloc(sizeof(double*), darp->N)));
	for (i=0; i < darp->N; i++)
		assert((darp->t[i] = calloc(sizeof(double), darp->N)));

	for (i=0; i < darp->N; i++) {
		int id;
		// <i> <x coord xi> <y coord yi> <service time at node di> <num passengers qi> <earliest time ei> <latest time li>
		err = scanf("%d %f %f %f %f %f %f\n", &id, &darp->x[i], &darp->y[i], &darp->d[i], &darp->q[i], &darp->e[i], &darp->l[i]);
		if (err != 7) {
			fprintf(stderr, "error processing line %d\n", i+2);
			exit(EXIT_FAILURE);
		}

		for (j=0; j <= i; j++)
			darp->t[i][j] = sqrt(pow(darp->x[j]-darp->x[i],2) + pow(darp->y[j]-darp->y[i],2));
	}
}


int
main(void)
{
	int restart, maxrestarts = 100;
	float cost_star = RAND_MAX;
	float *tardvec_star;
	float *s, *tardvec, tardmax, tardsum, ***arcmat;
	int i, j, *priority_rule;

	darp_t darp;
	aux_t aux;

	darp_create(&darp);
	initialize(&aux, darp);

	assert((s = calloc(sizeof(float), darp.N)));
	assert((priority_rule = calloc(sizeof(int), darp.N)));
	assert((tardvec = calloc(sizeof(float), darp.N)));
	assert((tardvec_star = calloc(sizeof(float), darp.N)));
	assert((arcmat = calloc(sizeof(int*), darp.N)));
	for (i=0; i < darp.N; i++) {
		assert((arcmat[i] = calloc(sizeof(int*), darp.N)));
		for (j=0; j < darp.N; j++)
			assert((arcmat[i][j] = calloc(sizeof(int), darp.m)));
	}

	for (restart=0; restart < maxrestarts; restart++) {
		prioritize(&aux, darp, tardvec, priority_rule);
		rho(&aux, darp, priority_rule, s, arcmat, tardvec, &tardmax, &tardsum);

		if (tardmax < cost_star) {
			printf("cost %f tardmax %f\n", tardsum, tardmax);
			cost_star = tardmax;
			memcpy(tardvec_star, tardvec, sizeof(float)*darp.N);
		} else {
			memcpy(tardvec, tardvec_star, sizeof(float)*darp.N);
		}
	}

	return 0;
}
