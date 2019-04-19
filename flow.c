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
	int *a, *o;
	float *s;
	int *v;
	float *tardiness, *best_objective_tardiness;
	int *inserted;
	int ***arc;
	tuple_t *tuples;
	chain_node_t **chain;
	float objective, max_tard;
} work_t;

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
work_create(work_t *work, darp_t darp)
{
	int i, j, k;

	work->prioritize = priority_criterion;

	assert((work->tuples = calloc(sizeof(tuple_t), darp.N)));
	assert((work->tardiness = calloc(sizeof(float), darp.N)));
	assert((work->s = calloc(sizeof(float), darp.N)));
	assert((work->v = calloc(sizeof(int), darp.n)));
	assert((work->best_objective_tardiness = calloc(sizeof(float), darp.N)));
	assert((work->inserted = calloc(sizeof(int), darp.N)));
	assert((work->arc = calloc(sizeof(int*), darp.N)));
	for (i=0; i < darp.N; i++) {
		assert((work->arc[i] = calloc(sizeof(int*), darp.N)));
		for (j=0; j < darp.N; j++)
			assert((work->arc[i][j] = calloc(sizeof(int), darp.m)));
	}
	assert((work->a = calloc(sizeof(int), darp.N)));
	assert((work->o = calloc(sizeof(int), darp.N)));
	assert((work->chain = calloc(sizeof(chain_node_t*), darp.m)));
	for (k=0; k < darp.m; k++) 
		assert((work->chain[k] = calloc(sizeof(chain_node_t), darp.N)));

}

static void
work_prioritize(work_t *work, darp_t darp)
{
	int i, z;
	// sort according to priority criterion
	for (i=0; i < darp.N; i++) {
		work->tuples[i].i = i;
		work->tuples[i].e = darp.e[i];
		work->tuples[i].l = darp.l[i] - 0.5*work->tardiness[i]*rand()/RAND_MAX;
		work->tuples[i].source = (i == 0);
		work->tuples[i].sink = (i == darp.N-1);
	}
	qsort(work->tuples, darp.N, sizeof(tuple_t), work->prioritize);

	// initialize priority rule
	for (z=0; z < darp.N; z++) {
		work->a[z] = work->tuples[z].i;
		work->o[work->tuples[z].i] = z;
	}
}

static void
work_schedule(work_t *work, darp_t darp)
{
	int i, j, k, p, kmin, kmax;
	int num_inserted;
	float tard, best_tard;
	int best_k, best_found;

	// initialize chain-form schedule
	for (k=0; k < darp.m; k++) {
		work->chain[k][0].i = 0;
		work->chain[k][0].leak = darp.Q;
	}

	// generate schedule
	work->s[0] = 0;
	num_inserted = 0;
	memset(work->inserted, 0, sizeof(int)*darp.N);
	memset(work->tardiness, 0, sizeof(float)*darp.N);

	work->objective = 0;
	work->max_tard = 0;
	while (num_inserted < darp.N-2)
	{
		// pick next eligible event from priority rule a[]
		for (p=0; p < darp.N; p++) {
			j = work->a[p];

			if (j == 0 || j == darp.N-1) 
				continue;

			//printf("\n%d ", j);

			if (work->inserted[j]) {
				//printf("already inserted");
				continue;
			}

			if (j > darp.n && j < 2*darp.n+1) {
				if (!work->inserted[j-darp.n]) {
					//printf("pickup not yet inserted");
					continue;
				}
				kmin = kmax = work->v[j-darp.n];
			} else {
				kmin = 0;
				kmax = darp.m-1;
			}

			best_tard = 10*darp.T;
			best_found = best_k = 0;

			for (k=kmin; k <= kmax; k++) {
				if (work->chain[k][0].leak < darp.q[j]) 
					continue;

				i = work->chain[k][0].i;
				work->s[j] = work->s[i] + darp.t[i][j] + darp.d[i];
				work->s[j] = fmax(darp.e[j], work->s[j]);
				tard = work->s[j] - darp.l[j];

				if (tard <= best_tard) {
					best_k = k;
					best_tard = tard;
					best_found = 1;
				}
			}

			if (best_found) {
				work->objective += darp.t[work->chain[best_k][0].i][j];
				work->max_tard = fmax(work->max_tard, best_tard);
				//printf("(%d,%d):%d\n", chain[best_k][0].i, j, k);
				work->arc[work->chain[best_k][0].i][j][best_k] = 1;
				work->chain[best_k][0].i = j;
				work->chain[best_k][0].leak -= darp.q[j];
				work->inserted[j] = 1;
				work->tardiness[j] = best_tard;
				if (j <= darp.n)
					work->v[j] = best_k;
				num_inserted++;
				//printf("inserted %d (k=%d leak=%d tard=%.1f)\n", j, best_k, chain[best_k][0].leak, best_tard);
			} else {
				//printf("no sufficient leak");
			}
		}
	}

	// append sink to all machines
	for (k=0; k < darp.m; k++) {
		i = work->chain[k][0].i;
		work->arc[i][darp.N-1][k] = 1;
		work->s[darp.N-1] = fmax(work->s[darp.N-1], work->s[i] + darp.d[i] + darp.t[i][darp.N-1]);
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
	int restart, maxrestarts = 200;
	float best_objective = RAND_MAX;
	float *best_objective_tardiness;

	darp_t darp;
	work_t work;

	darp_create(&darp);
	work_create(&work, darp);

	assert((best_objective_tardiness = calloc(sizeof(float), darp.N)));

	for (restart=0; restart < maxrestarts; restart++) {

		work_prioritize(&work, darp);
		work_schedule(&work, darp);

		if (work.max_tard < best_objective) {
			printf("objective %f max_tard %f\n", work.objective, work.max_tard);
			best_objective = work.max_tard;
			memcpy(best_objective_tardiness, work.tardiness, sizeof(float)*darp.N);
		} else {
			memcpy(work.tardiness, best_objective_tardiness, sizeof(float)*darp.N);
		}
	}

	return 0;
}
