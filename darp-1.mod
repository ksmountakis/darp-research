param n integer >= 1;
param m integer >= 1;
param T integer >= 0;
param L integer >= 0;
param Q integer >= 0;

set N := 1 .. 2*n+2;
set P := 2 .. n+1;
set D := n+2 .. 2*n+1;
set PD := 2 .. 2*n+1;
set K := 1 .. m;

param d {i in N} integer >= 0;
param q {i in N} integer;
param e {i in N} >= 0;
param l {i in N} >= 0;
param t {i in N, j in N} >= 0;
param bigM integer := 1e3*T;

var x {i in N, j in N, k in K} binary;
#var s {i in N, k in K} >= 0;
var s {i in N} >= 0;
var u {i in N} >= 0;

var f;

minimize cost: f;

subject to objective:
	f = sum {i in N, j in N, k in K} t[i,j]*x[i,j,k];
	#f = 100*((sum {i in P} ((s[i+n] - s[i])/t[i,i+n]))/n - 1.0);

subject to outP{i in P}:
	sum {j in N, k in K} x[i,j,k] = 1;

subject to outP_outD{i in P, k in K}:
	sum {j in N} x[i,j,k] = sum {j in N} x[i+n, j, k];

subject to in_and_out {i in PD, k in K}:
	sum {j in N} x[j,i,k] = sum {j in N} x[i,j,k];

subject to source {k in K}:
	sum {j in P} x[1,j,k] = 1;

subject to source_no_D {k in K}:
	sum {j in D} x[1,j,k] = 0;

subject to sink {k in K}:
	sum {i in D} x[i,2*n+2,k] = 1;

subject to sink_no_P {k in K}:
	sum {i in P} x[i,2*n+2,k] = 0;


subject to precedence_no_cycles {i in N, j in N}:
	s[j] >= s[i] + d[i] + t[i,j] - (1 - sum {k in K} x[i,j,k])*bigM;
#subject to precedence_no_cycles {i in N, j in N, k in K}:
#	s[j,k] >= s[i,k] + d[i] + t[i,j] - (1 - x[i,j,k])*bigM;
#


subject to resource_use {i in N, j in N}:
	u[j] >= u[i] + q[i] - (1 - sum {k in K} x[i,j,k])*bigM;

subject to capacity {i in N}:
	u[i] <= Q;

subject to earliest_start {i in N}:
	s[i] >= e[i];
#subject to earliest_start {i in N, k in K}:
#	s[i, k] >= e[i];

subject to latest_start {i in N}:
	s[i] <= l[i];
#subject to latest_start {i in N, k in K}:
#	s[i,k] <= l[i];

subject to maxusertime {i in P}:
	s[i+n] - (s[i] + d[i]) <= L;
#subject to maxusertime {i in P, k in K}:
#	s[i+n, k] - (s[i,k] + d[i]) <= L;

subject to minusertime {i in P}:
	s[i+n] - (s[i] + d[i]) >= t[i,i+n];
#subject to minusertime {i in P, k in K}:
#	s[i+n, k] - (s[i, k] + d[i]) >= t[i,i+n];

subject to maxridetime:
	s[2*n+2] - s[1] <= T;
#subject to maxridetime {k in K}:
#	s[2*n+2, k] - s[1, k] <= T;

solve;
display f, u;

