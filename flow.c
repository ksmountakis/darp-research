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

int m, n, N;

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

int
main(void)
{
	int err;
	float T, Q, L;
	int i, j, k, z;
	float *x, *y, *d, *q, *e, *l;
	double **t;
	int *a, *o;

	int num_inserted;
	float *s;
	int *v;
	float *tardiness, *best_objective_tardiness;
	int *inserted;
	int ***arc;
	int p;
	float tard, best_tard;
	int best_i, best_k, best_found;
	float objective;
	float max_tard;
	int kmin, kmax;
	tuple_t *tuples;
	chain_node_t **chain;
	int num_restarts, max_restarts = 1000;
	float best_objective = RAND_MAX;

	// <m> <n> <max veh ride minutes Tk> <veh capacity Qk> <max user ride minutes Li>
	err = scanf("%d %d %f %f %f\n", &m, &N, &T, &Q, &L);
	if (err != 5) {
		fprintf(stderr, "error processing line 1\n");
		exit(EXIT_FAILURE);
	}

	N = 2*N + 2;
	n = (N-2)/2;

	fprintf(stderr, "n %d m %d\n", n, m);

	// allocations
	assert((x = calloc(sizeof(float), N)));
	assert((y = calloc(sizeof(float), N)));
	assert((d = calloc(sizeof(float), N)));
	assert((q = calloc(sizeof(float), N)));
	assert((e = calloc(sizeof(float), N)));
	assert((l = calloc(sizeof(float), N)));
	assert((t = calloc(sizeof(double*), N)));
	for (i=0; i < N; i++)
		assert((t[i] = calloc(sizeof(double), N)));

	assert((s = calloc(sizeof(float), N)));
	assert((v = calloc(sizeof(int), n)));
	assert((tardiness = calloc(sizeof(float), N)));
	assert((best_objective_tardiness = calloc(sizeof(float), N)));
	assert((inserted = calloc(sizeof(int), N)));
	assert((arc = calloc(sizeof(int*), N)));
	for (i=0; i < N; i++) {
		assert((arc[i] = calloc(sizeof(int*), N)));
		for (j=0; j < N; j++)
			assert((arc[i][j] = calloc(sizeof(int), m)));
	}
	assert((tuples = calloc(sizeof(tuple_t), N)));
	assert((a = calloc(sizeof(int), N)));
	assert((o = calloc(sizeof(int), N)));
	assert((chain = calloc(sizeof(chain_node_t*), m)));
	for (k=0; k < m; k++) 
		assert((chain[k] = calloc(sizeof(chain_node_t), N)));
	
	for (i=0; i < N; i++) {
		int id;
		// <i> <x coord xi> <y coord yi> <service time at node di> <num passengers qi> <earliest time ei> <latest time li>
		err = scanf("%d %f %f %f %f %f %f\n", &id, &x[i], &y[i], &d[i], &q[i], &e[i], &l[i]);
		if (err != 7) {
			fprintf(stderr, "error processing line %d\n", i+2);
			exit(EXIT_FAILURE);
		}

		for (j=0; j <= i; j++)
			t[i][j] = sqrt(pow(x[j]-x[i],2) + pow(y[j]-y[i],2));
	}

	num_restarts = 0;
restart:
	// sort according to priority criterion
	for (i=0; i < N; i++) {
		tuples[i].i = i;
		tuples[i].e = e[i];
		//tuples[i].l = l[i] - (0.5*rand()/RAND_MAX*tardiness[i]);
		tuples[i].l = l[i] - 0.5*tardiness[i]*rand()/RAND_MAX;
		tuples[i].source = (i == 0);
		tuples[i].sink = (i == N-1);
	}
	qsort(tuples, N, sizeof(tuple_t), priority_criterion);

	// initialize priority rule
	//printf("a: ");
	for (z=0; z < N; z++) {
		a[z] = tuples[z].i;
		o[tuples[z].i] = z;
		//printf("%2d ", a[z]);
	}
	//puts("");

	// initialize chain-form schedule
	for (k=0; k < m; k++) {
		chain[k][0].i = 0;
		chain[k][0].leak = Q;
	}

	// generate schedule
	s[0] = 0;
	num_inserted = 0;
	memset(inserted, 0, sizeof(int)*N);
	memset(tardiness, 0, sizeof(float)*N);
	objective = 0;
	max_tard = 0;

	while (num_inserted < N-2)
	{
		// pick next eligible event from priority rule a[]
		for (p=0; p < N; p++) {
			j = a[p];

			if (j == 0 || j == N-1) 
				continue;

			//printf("\n%d ", j);

			if (inserted[j]) {
				//printf("already inserted");
				continue;
			}

			if (j > n && j < 2*n+1) {
				if (!inserted[j-n]) {
					//printf("pickup not yet inserted");
					continue;
				}
				kmin = kmax = v[j-n];
			} else {
				kmin = 0;
				kmax = m-1;
			}

			best_tard = 10*T;
			best_i = best_found = best_k = 0;

			for (k=kmin; k <= kmax; k++) {
				if (chain[k][0].leak < q[j]) 
					continue;

				i = chain[k][0].i;
				s[j] = s[i] + t[i][j] + d[i];
				s[j] = fmax(e[j], s[j]);
				tard = s[j] - l[j];

				if (tard <= best_tard) {
					best_i = i;
					best_k = k;
					best_tard = tard;
					best_found = 1;
				}
			}

			if (best_found) {
				objective += t[chain[best_k][0].i][j];
				max_tard = fmax(max_tard, best_tard);
				//printf("(%d,%d):%d\n", chain[best_k][0].i, j, k);
				arc[chain[best_k][0].i][j][best_k] = 1;
				chain[best_k][0].i = j;
				chain[best_k][0].leak -= q[j];
				inserted[j] = 1;
				tardiness[j] = best_tard;
				if (j <= n)
					v[j] = best_k;
				num_inserted++;
				//printf("inserted %d (k=%d leak=%d tard=%.1f)\n", j, best_k, chain[best_k][0].leak, best_tard);
			} else {
				//printf("no sufficient leak");
			}
		}
	}

	// append sink to all machines
	for (k=0; k < m; k++) {
		i = chain[k][0].i;
		arc[i][N-1][k] = 1;
		s[N-1] = fmax(s[N-1], s[i]+d[i]+t[i][N-1]);
	}

#if 1
	// verify constraints
	for (i=0; i < N; i++) {
		int sum = 0;
		for (k=0; k < m; k++) {
			for (j=0; j < N; j++) {
				sum += arc[i][j][k];
			}
		}
		//printf("sum {j in N, k in K} arc[%d,j,k] (%d)\n", i, sum);
	}
#endif

	if (max_tard < best_objective) {
		printf("objective %f max_tard %f\n", objective, max_tard);
		best_objective = max_tard;
		memcpy(best_objective_tardiness, tardiness, sizeof(float)*N);
	} else {
		memcpy(tardiness, best_objective_tardiness, sizeof(float)*N);
	}

	if (num_restarts++ < max_restarts)
		goto restart;

#if 0
	// process each node in order
	printf("P: ");
	for (i=1; i <= n; i++) 
		printf("%2d:%2d ", i, o[i]);
	puts("");

	printf("D: ");
	for (i=1+n; i <= 2*n; i++) 
		printf("%2d:%2d ", i, o[i]);

	printf("%2d:%2d\n", N-1, o[N-1]);
#endif
	
	return 0;
}
