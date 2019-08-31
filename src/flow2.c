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
	int i;
	int leak;
	float s;
} ChainNode;

typedef struct {
	int N, n, m;
	float T, Q, L;
	int i, j, k, z;
	float *x, *y, *d, *q, *e, *l;
	float **t;
	ChainNode **chain;
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

static int
inD(Darp *darp, int i) {
	return (i > darp->n);
}

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
	for (i=1; i < darp->n; i++) {
		darp->l[i] = fmin(darp->l[i], darp->l[i+darp->n]-darp->t[i][i+darp->n]);
		if (inD(darp, i))
			darp->e[i]=0;
	}

	// allocate chains
	int k;
	assert((darp->chain = calloc(sizeof(ChainNode*),darp->m)));
	for (k=0; k < darp->m; k++) {
		assert((darp->chain[k] = calloc(sizeof(ChainNode),darp->N)));
		darp->chain[k][0].i = 0;
		darp->chain[k][0].leak = darp->Q;
	}
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

static float
darprho(Darp *darp, ChainNode **chain, RuleItem *a, int N, float *tardvec, float *tardmaxptr)
{
	int k;
	// initialize schedule
	a[0].s = 0;
	for (k=0; k < darp->m; k++) {
		chain[k][0].i = a[0].i;
		chain[k][0].s = a[0].s;
		chain[k][0].leak = darp->Q;
	}

	memset(tardvec, 0, sizeof(float)*darp->N);
	float tardmax = 0;
	float cost = 0;

	// mapping from node to {0,1}
	static int *inserted;
	assert((inserted = realloc(inserted, sizeof(int)*darp->N)));
	memset(inserted, 0, sizeof(int)*darp->N);
	inserted[0] = 1;
	int num_inserted = 1;

	// mapping from node to {0,...,m-1}
	static int *vehicle;
	assert((vehicle = realloc(vehicle, sizeof(int)*darp->N)));
	memset(vehicle, 0, sizeof(int)*darp->N);

	while (num_inserted < N)
	{
		int p;
		// pick next eligible event from priority rule a[]
		for (p=0; p < N; p++) 
		{
			if (inserted[a[p].i])
				continue;

			// determine range of vehicles to try
			int kmin, kmax;
			if (a[p].i > darp->n && a[p].i < 2*darp->n+1) {
				if (!inserted[a[p].i-darp->n]) {
					//printf("pickup not yet inserted");
					continue;
				}
				kmin = kmax = vehicle[a[p].i-darp->n];
			} else {
				kmin = 0;
				kmax = darp->m-1;
			}

			// try each vehicle
			int best_found = 0, best_k = 0;
			float best_tard = 1e9, best_s;

			for (k=kmin; k <= kmax; k++) 
			{
				float tard;
				if (chain[k][0].leak < darp->q[a[p].i]) 
					continue;

				a[p].s = fmax(darp->e[a[p].i], chain[k][0].s + darp->t[chain[k][0].i][a[p].i] + darp->d[a[p].i]);
				tard = a[p].s - darp->l[a[p].i];

				if (tard < best_tard) {
					best_k = k;
					best_tard = tard;
					best_found = 1;
					best_s = a[p].s;
				}
			}

			// finalize best position found
			if (best_found) {
				a[p].s = best_s;
				tardmax = fmax(tardmax, best_tard);
				tardvec[a[p].i] = best_tard;
				cost += darp->t[chain[best_k][0].i][a[p].i];
				//arcmat[rho_ctx->chain[best_k][0].i][j][best_k] = 1;
				//printf("+ (%2d -> %2d) k %d s %.1f e %.1f leak %.1f\n", chain[best_k][0].i, a[p].i, best_k, a[p].s, darp->e[a[p].i], chain[best_k][0].leak - darp->q[a[p].i]);
				chain[best_k][0].i = a[p].i;
				chain[best_k][0].leak -= darp->q[a[p].i];
				chain[best_k][0].s = a[p].s;
				inserted[a[p].i] = 1;
				if (a[p].i <= darp->n)
					vehicle[a[p].i] = best_k;
				num_inserted++;
			} else {
				//printf("position not found for: %2d\n", a[p].i);
			}
		}
	}

	*tardmaxptr = tardmax;

	return cost;
}

static float
darpinsert(Darp *darp, Schedule *s, int j)
{
	int lmin, l, y, i, accept;
	RuleItem a[darp->N];
	static float *tardvec;
	float tardmax, cost;

	assert((tardvec = realloc(tardvec, sizeof(float)*darp->N)));

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
	if (j == darp->N-1)
		lmin = s->n;

	printf("inserting %c %2d: ", PDchar(darp,j), j);
	// print current schedule
	for (l=0; l < s->n; l++) 
		printf("%2d ", s->a[l].i);
	printf(" lmin=%2d\n", lmin);


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
		cost = darprho(darp, darp->chain, a, s->n+1, tardvec, &tardmax);

		// examine schedule
		accept = 0;
		if (tardmax <= 0.0 || l == s->n)
			accept = 1;

		printf("a: ");
		for (i=0; i <= s->n; i++) {
			float dt = a[i].s - darp->l[a[i].i];
			printf("%2d ", a[i].i);
			if (dt > 0)
				printf("[%2.1f] ", dt);
			else
				printf(" ");
		}
		printf(" := %.1f\n", cost);

		if (accept) 
			break;
	}


	assert(accept);
	s->n++;
	for (y=0; y < s->n; y++) 
		s->a[y] = a[y];

#if 0
	fprintf(stderr, "a(+%2d): ", j);
	for (i=0; i <= s->n; i++) {
		float dt = s->a[i].s - darp->l[s->a[i].i];
		fprintf(stderr, "%2d ", s->a[i].i);
		if (dt > 0)
			fprintf(stderr, "[%2.1f] ", dt);
		else
			fprintf(stderr, " ");
	}
	fprintf(stderr, " := %.1f\n", cost);
#endif

	return cost;
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

		printf("RULE: ");
		for (i=0; i < darp.N; i++)
			printf("%2d ", rule[i].i);
		puts("");

		// create schedule
		darpreset(&darp, &schedule);
		for (i=1; i < darp.N; i++) {
			if (inD(&darp, rule[i].i))
				continue;
			darpinsert(&darp, &schedule, rule[i].i);
			darpinsert(&darp, &schedule, rule[i].i+darp.n);
		}

		// examine schedule
		// TODO
	}

	return 0;
}
