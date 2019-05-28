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

param minstep := 1.0;
param t_int {i in N, j in N} := round(t[i,j] / minstep);
param bigM integer := ceil(T/minstep);
param TMAX := ceil(T/minstep);
set Z := 0 .. TMAX;

var x {i in N, z in Z, k in K} binary;
var X {i in N, z in Z} >= 0;
var s {i in N} >= 0;
 
var f;

minimize cost: f;

subject to objective:
	f = sum {i in P} s[i];

# define X[i,z] 
subject to define_Xiz {i in N, z in Z}:
	X[i,z] = sum {k in K} x[i,z,k];

# define s[i]
subject to define_Si {i in N}:
	s[i] = sum {z in Z} z * X[i,z];

# at most one event per time-point per time-line
subject to atmost_one_event_per {z in Z, k in K}:
	sum {i in N} x[i,z,k] <= 1;

# exactly one time-point per event
subject to exactly_one_timepoint_per {i in N}:
	sum {z in Z} X[i, z] = 1;

# drop-offs on same time-line as pick-ups
subject to same_timeline_per {i in D, k in K}:
	sum {z in Z} x[i,z,k] = sum {z in Z} x[i-n,z,k];

# drop-offs after pick-ups
subject to dropoff_after_pickup_per {i in P}:
	s[i+n] >= s[i] + t_int[i,i+n];

# pair-wise distance constraints
subject to distance_constraint_per {i in N, j in N, z0 in 1 .. TMAX}:
	s[j] >= s[i] + t_int[i,j] + d[i] + (X[j,z0] - sum {z in 0 .. z0-1} X[i,z] - 2)*bigM;

# vehicle capacity constraints
subject to capacity_constraint_per {k in K, z0 in Z}:
	sum {z in 0 .. z0, j in N} q[j]*x[j,z,k] <= Q;

# max user time
subject to maxusertime {i in P}:
	s[i+n] - (s[i] + round(d[i]/minstep)) <= ceil(L/minstep);

# max vehicle time
subject to maxvehtime {k in K}:
	sum {z in Z, i in N} z * x[i,z,k] <= ceil(T/minstep);

# earliest time
subject to earliest_start {i in N}:
	s[i] >= floor(e[i]/minstep);

# latest time
subject to latest_start {i in N}:
	s[i] <= ceil(l[i]/minstep);

subject to s0: 
	s[1] = 0;

solve;
display s;
