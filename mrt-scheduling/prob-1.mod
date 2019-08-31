param n integer >= 1;
param d {1 .. n} >= 0;
param lower {1 .. n} >= 0 default 0;
param upper {1 .. n} >= 0 default 10e5;
param mrt {1 .. n, 1 .. n} >= 0 default 10e5;

var t {1 .. n} >= 0;

minimize cost:
    sum {i in 1 .. n} t[i];

subject to sequencing {i in 1 .. n-1}:
    t[i] - t[i+1] <= -d[i];

subject to lower_bounds {i in 1 .. n}:
    t[i] >= lower[i];

subject to upper_bounds {i in 1 .. n}:
    t[i] <= upper[i];

subject to max_ride_bounds {i in 1 .. n, j in i .. n}:
    t[j] - t[i] <= mrt[i,j];


#option solver kestrel;
#option kestrel_options "solver=CPLEX";

solve;

display cost, t;
