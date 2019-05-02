#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#if 0

typedef struct {
	int i;
	int leak;
	float s;
} chain_node_t;

typedef struct {
	int (*make_priority_rule)(const void*, const void*);
	int *o;
	float *s;
	int *v;
	int *inserted;
	tuple_t *tuples;
	chain_node_t **chain;
} rho_ctx_t;
#endif

typedef struct {
	int N, n, m;
	float T, Q, L;
	int i, j, k, z;
	float *x, *y, *d, *q, *e, *l;
	float **t;
} Darp;

typedef struct {
	int i;
	float s;
	int source;
	int sink;
} RuleItem;

typedef struct {
	int n;
	RuleItem *a;
} Schedule;

static void
darpread(Darp *darp)
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
	assert((darp->t = calloc(sizeof(float*), darp->N)));
	for (i=0; i < darp->N; i++)
		assert((darp->t[i] = calloc(sizeof(float), darp->N)));

	// read x,y locations
	for (i=0; i < darp->N; i++) {
		int id;
		// <i> <x coord xi> <y coord yi> <service time at node di> <num passengers qi> <earliest time ei> <latest time li>
		err = scanf("%d %f %f %f %f %f %f\n", &id, &darp->x[i], &darp->y[i], &darp->d[i], &darp->q[i], &darp->e[i], &darp->l[i]);
		if (err != 7) {
			fprintf(stderr, "error processing line %d\n", i+2);
			exit(EXIT_FAILURE);
		}
	}

	// derive pair-wise distances
	for (i=0; i < darp->N; i++) 
		for (j=0; j < darp->N; j++) 
			darp->t[i][j] = sqrtf(powf(darp->x[j]-darp->x[i],2) + powf(darp->y[j]-darp->y[i],2));

	// tighten windows
	for (i=1; i < darp->n; i++) 
		darp->l[i] = fmin(darp->l[i], darp->l[i+darp->n]-darp->t[i][i+darp->n]);
}

static int 
darprulecmp(const void *t1ptr, const void *t2ptr) {
	RuleItem *t1 = (RuleItem*)t1ptr;
	RuleItem *t2 = (RuleItem*)t2ptr;
	float time1 = t1->s;
	float time2 = t2->s;

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

static int
inD(Darp *darp, int i) {
	return (i > darp->n);
}

static int*
darpmkrule(Darp *darp, RuleItem **rule)
{
	static int *order;

	assert((*rule = realloc(*rule, sizeof(RuleItem)*(darp->N))));
	assert((order = realloc(order, sizeof(int)*(darp->N))));

	int i;
	// sort 
	for (i=0; i < darp->N; i++) {
		(*rule)[i].i = i;
		(*rule)[i].s = darp->l[i];
		(*rule)[i].source = (i == 0);
		(*rule)[i].sink = (i == darp->N-1);
	}
	qsort((*rule), darp->N, sizeof(RuleItem), darprulecmp);

	for (i=0; i < darp->N; i++) 
		order[(*rule)[i].i] = i;

	// sanity check
	for (i=1; i <= darp->n; i++)
		if (order[i] > order[i+darp->n])
			fprintf(stderr, "darpmkrule: D %d precedeces P %d\n", i+darp->n, i);

	return order;
}


#if 0
static void
rho_ctx_create(darp_t darp, rho_ctx_t *rho_ctx)
{
	int k;

	rho_ctx->make_priority_rule = priority_criterion;
	assert((rho_ctx->tuples = calloc(sizeof(tuple_t), darp.N)));
	assert((rho_ctx->v = calloc(sizeof(int), darp.n)));
	assert((rho_ctx->inserted = calloc(sizeof(int), darp.N)));
	assert((rho_ctx->o = calloc(sizeof(int), darp.N)));
	assert((rho_ctx->chain = calloc(sizeof(chain_node_t*), darp.m)));
	for (k=0; k < darp.m; k++) 
		assert((rho_ctx->chain[k] = calloc(sizeof(chain_node_t), darp.N)));
}

static void
rho(rho_ctx_t *rho_ctx, const darp_t darp, int *a, float *s, float ***arcmat, float *tardvec, float *tardmax, float *tardsum)
{
	int i, j, k, p, kmin, kmax;
	int num_inserted;
	float tard, best_tard;
	int best_k, best_found;

	// initialize chain-form schedule
	for (k=0; k < darp.m; k++) {
		rho_ctx->chain[k][0].i = 0;
		rho_ctx->chain[k][0].leak = darp.Q;
	}

	// generate schedule
	s[0] = 0;
	num_inserted = 0;
	memset(rho_ctx->inserted, 0, sizeof(int)*darp.N);
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

			if (rho_ctx->inserted[j]) {
				//printf("already inserted");
				continue;
			}

			if (j > darp.n && j < 2*darp.n+1) {
				if (!rho_ctx->inserted[j-darp.n]) {
					//printf("pickup not yet inserted");
					continue;
				}
				kmin = kmax = rho_ctx->v[j-darp.n];
			} else {
				kmin = 0;
				kmax = darp.m-1;
			}

			best_tard = 10*darp.T;
			best_found = best_k = 0;

			for (k=kmin; k <= kmax; k++) {
				if (rho_ctx->chain[k][0].leak < darp.q[j]) 
					continue;

				i = rho_ctx->chain[k][0].i;
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
				*tardsum += darp.t[rho_ctx->chain[best_k][0].i][j];
				*tardmax = fmax(*tardmax, best_tard);
				//printf("(%d,%d):%d\n", chain[best_k][0].i, j, k);
				arcmat[rho_ctx->chain[best_k][0].i][j][best_k] = 1;
				rho_ctx->chain[best_k][0].i = j;
				rho_ctx->chain[best_k][0].leak -= darp.q[j];
				rho_ctx->inserted[j] = 1;
				tardvec[j] = best_tard;
				if (j <= darp.n)
					rho_ctx->v[j] = best_k;
				num_inserted++;
				//printf("inserted %d (k=%d leak=%d tard=%.1f)\n", j, best_k, chain[best_k][0].leak, best_tard);
			} else {
				//printf("no sufficient leak");
			}
		}
	}

	// append sink to all machines
	for (k=0; k < darp.m; k++) {
		i = rho_ctx->chain[k][0].i;
		arcmat[i][darp.N-1][k] = 1;
		s[darp.N-1] = fmax(s[darp.N-1], s[i] + darp.d[i] + darp.t[i][darp.N-1]);
	}
}
#endif

#define PDchar(darp, j) ((inD(darp,j))? 'D': 'P')

static void
darpreset(Darp *darp, Schedule *s)
{
	assert((s->a = realloc(s->a, sizeof(RuleItem)*darp->N)));
	memset(s->a, 0, sizeof(RuleItem)*darp->N);
	s->n = 1;
	s->a[0].i = 0;
	s->a[s->a[0].i].s = 0;
	s->a[s->a[0].i].source = 1;
}


static void
darpinsert(Darp *darp, Schedule *s, int j)
{
	int lmin, l, y, i, accept;
	RuleItem a[darp->N];

	// compute s->a[] by sorting in ascending start-time
	qsort(s->a, s->n, sizeof(RuleItem), darprulecmp);

	// compute lmin
	lmin = 1;
	if (inD(darp, j)) {
		lmin = -1;
		for (l=0; l < s->n; l++) {
			if (s->a[l].i == j-darp->n) {
				lmin = l+1;
				break;
			}
		}
		if(lmin == -1) {
			fprintf(stderr, "darpinsert(%d): P %d not already inserted\n", j, j-darp->n);
			exit(1);
		}
	}

	printf("inserting %c %2d: ", PDchar(darp,j), j);
	for (l=0; l < s->n; l++)
		printf("%2d ", s->a[l].i);
	printf(" lmin=%2d\n", lmin);

	if (j == darp->N-1)
		lmin = s->n;

	memset(a, 0, sizeof(RuleItem)*darp->N);
	for (l=lmin; l <= s->n; l++) {
		// prepare a[] with j as the l-th element
		for (y=0; y <= l-1; y++) 
			a[y] = s->a[y];

		a[l].i = j;
		a[l].s = -1;
		a[l].source = (j == 0);
		a[l].sink = (j == darp->N-1);

		for (y=l; y <= s->n-1; y++) 
			a[y+1] = s->a[y];

		// create s[] based on a[]
		// TODO
		a[l].s = (float)j;
		
		// examine schedule
		accept = 0;
		if (l == s->n)
			accept = 1;

		printf("a: ");
		for (i=0; i <= s->n; i++)
			printf("%2d ", a[i].i);
		printf("\n");

		if (accept) 
			break;
	}

	assert(accept);
	s->n++;
	for (y=0; y < s->n; y++) 
		s->a[y] = a[y];
}

int
main(void)
{
	int i, loop, maxloop = 1;
	Darp darp;
	RuleItem *rule = NULL;
	Schedule schedule = {0};

	// read problem
	darpread(&darp);

	// perform main iterations
	for (loop=0; loop < maxloop; loop++) {
		// create priority rule
		darpmkrule(&darp, &rule);

		// create schedule
		// TODO: reset schedule
		darpreset(&darp, &schedule);
		for (i=1; i < darp.N; i++) 
			darpinsert(&darp, &schedule, rule[i].i);

		// examine schedule
		// TODO
	}

	return 0;
}
