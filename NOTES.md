VRP Research Notes



Write the GAMA out:

Objective function:
f(X)=Cx where C and X are integer vectors
    flatten X and C matrices
    X contains both x and y variables
    we only need nlogn variables to represent y, we store bits of y as separate variables

Feasible Function Generator:

    A matrix: I might just have to code this one out it'll probably be hard to 


    It looks like I can just reuse sim_anneal_paths.py for the actual annealing for feasible functions, I just have to figure out how to find A, b, E, and L


CVRP Constrints:
Cost matrix is nonnegative and triangle inequality holds (c_ij<=c_ik+c_kj), c_ij=c_ji (more of a generation constraint)
All vehicles start and end at depot, so (i=0 and i=n+1 represent depot)
Notably, each vertex can be visited only once across all vehicles
Note:d_i<=q

Rewrite 14
y_j>=y_i+d_j*x_ij-Q(1-x_ij)
y_j-y_i-d_j*x_ij+Q(1-x_ij)>=0
y_j-y_i-d_j*x_ij+Q-Qx_ij>=0
y_j-y_i-(d_j+Q)*x_ij+Q>=0
y_j-y_i-(d_j+Q)*x_ij+Q-S=0
y_j-y_i-(d_j+Q)*x_ij-S=-Q
LHS i think is bounded by 2Q so we need n*(n-1)(log(Q)+1)