from generation.generate_cost_matrix import generate_cost_matrix
from generation.generate_graph import generate_graph
from compute_feasible_sol_matrices import compute_matrices
from sim_anneal_paths import get_feasible

import numpy as np
import time
import subprocess


#Configure this:
#n is # of nodes in graph
n=10
#q_clamp is a constant multiplied to each random number of d[i], higher q_clamp results in a graph
#with more feasible solutions (since there are more possible paths each car can take before it runs
#out of capacity)
q_clamp=1
#SA_samples can be increased if not enough feasible solutions are being found
SA_samples=10_000
#max_gbasis_size dictates how large the partial graver basis can be, can be increased for better answers
#but 100k is already quite good
max_gbasis_size=100_000
#only used if there are more feasible solutions than gama_num_seeds, picks number of random feasible solutions to
#start with
gama_num_seeds=1_000



start_time=time.time()
cost = generate_cost_matrix(n)
np.savetxt('data/cost_matrix.txt', cost, fmt="%.4f")
cmatrix_time = time.time()
print(f"Generated {n}x{n} cost matrix in data/cost_matrix.txt, cumulated time: {cmatrix_time-start_time}")
q,d,k = generate_graph(n,q_clamp)
graphgen_time = time.time()
print(f"Generated graph, cumulated time: {graphgen_time-start_time}")
A,b=compute_matrices(k,d,q)
SAmatrix_time = time.time()
print(f"SA matrices computed, cumulated time: {SAmatrix_time-start_time}")
get_feasible(A,b,len(d),q,d,SA_samples)
SAfin_time = time.time()
print(f"feasible solutions computed, cumulated time: {SAfin_time-start_time}")
result = subprocess.run(['./compute_graver_basis',str(n),str(max_gbasis_size)], capture_output=True, text=True)
g_basis_shape = result.stdout.strip().splitlines()[-1].split(' ')
print("Final g_basis size: "+g_basis_shape[0]+" x "+g_basis_shape[1])
gbasis_time = time.time()
print(f"graver basis computed, cumulated time: {gbasis_time-start_time}")
result = subprocess.run(['./gama_solver',str(n),str(g_basis_shape[0]),str(g_basis_shape[1]),str(gama_num_seeds)])





