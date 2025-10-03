VRP Research Notes



Write the GAMA out:

Objective function:\n
f(X)=Cx where C and X are integer vectors
    flatten X and C matrices
    X contains both x and y variables
    we only need nlogn variables to represent y, we store bits of y as separate variables\n

It looks like I can just reuse sim_anneal_paths.py for the actual annealing for feasible functions, I just have to figure out how to find A and b\n


CVRP Constrints:
Cost matrix is nonnegative and triangle inequality holds (c_ij<=c_ik+c_kj), c_ij=c_ji (more of a generation constraint)
All vehicles start and end at depot, so (i=0 and i=n+1 represent depot)\n
Notably, each vertex can be visited only once across all vehicles\n
Note:d_i<=q\n

Rewrite 14\n
y_j>=y_i+d_j\*x_ij-Q(1-x_ij)\n
y_j-y_i-d_j\*x_ij+Q(1-x_ij)>=0\n
y_j-y_i-d_j\*x_ij+Q-Qx_ij>=0\n
y_j-y_i-(d_j+Q)\*x_ij+Q>=0\n
y_j-y_i-(d_j+Q)\*x_ij+Q-S=0\n
y_j-y_i-(d_j+Q)\*x_ij-S=-Q\n
LHS i think is bounded by 2Q so we need n\*(n-1)(log(Q)+1)