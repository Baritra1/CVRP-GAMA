#adapted from https://doi.org/10.1287/ijoc.2024.0574.cd
import numpy as np
import pandas as pd
import dimod
import manual_sa.sampler as samplers
import neal
#check for subtours or vehicle taking more than capacity
#this is implemented in python because computing the graver basis is more of a bottleneck than this but I 
#could easily port to cpp if I find this starts taking too long for larger n
def check_possible(feas_sol,n,x,cum_demand,q,d):
    sum_visited=0
    if cum_demand>q:
        return sum_visited+1,False
    for i in range(n+1):
        if(i==n):
            return sum_visited+1,True
        if(feas_sol[x*(n+1)+i]==1):
            new_visit,poss = check_possible(feas_sol,n,i,cum_demand+d[i],q,d)
            sum_visited+=new_visit
            if(poss==False):
               return sum_visited+1,False    
    return sum_visited+1,True
                        


def get_feasible(A, b, n, q,d,samples=10000):

    AA = np.dot(A.T, A)
    h = -2.0*np.dot(b.T, A)
    Q = AA + np.diag(h)
    offset = np.dot(b.T, b) + 0.0
    bqm_model = dimod.BinaryQuadraticModel.from_numpy_matrix(mat=Q, offset=offset)
    simAnnSampler = samplers.SimulatedAnnealingSampler()
    sampler = simAnnSampler
    response = sampler.sample(bqm_model, num_reads=samples)
    response = response.aggregate()
    np.set_printoptions(threshold=np.inf)
    with open("data/raw_sols.txt", "w") as f:
        f.write(str(response.record)) # Add a newline character for each item
    filter_idx = [i for i, e in enumerate(response.record.energy) if e == 0.0]
    feas_sols = response.record.sample[filter_idx]
    feas_sols_filter=[]
    print("# of SA Solutions: "+str(len(feas_sols)))
    for i in feas_sols:
        visited,poss=check_possible(i,n,0,0,q,d)
        if(poss and visited==n):
            feas_sols_filter.append(i)
    feas_sols_clean=[]
    for i in feas_sols_filter:
        feas_sols_clean.append(i[:n*n+n])
    feas_sols_uniq = np.unique(feas_sols_clean, axis=0)
    print("# of possible unique SA Solutions "+str(len(feas_sols_uniq)))
    np.savetxt('data/feas_sols_sorted_2.txt', feas_sols_uniq)



def get_feasible_encode(A, b, E, L, sampler, samples=20):

    Q_I = np.dot(A.T, A)
    h = 2.0*np.dot(np.dot(L.T, Q_I) - np.dot(b.T, A), E)
    Q = np.dot(E.T, np.dot(Q_I, E)) + np.diag(h[0])
    offset = np.dot(L.T, np.dot(Q_I, L)) + np.dot(b.T, b) - 2.0*np.dot(b.T, np.dot(A, L))

    # Define Binary Quadratic Model
    bqm = dimod.BinaryQuadraticModel.from_numpy_matrix(mat=Q, offset=offset)

    response = sampler.sample(bqm, num_reads=samples)

    response = response.aggregate()

    filter_idx = [i for i, e in enumerate(response.record.energy) if e == 0.0]

    feas_sols_X = response.record.sample[filter_idx]

    print(L.shape, E.shape, feas_sols_X.shape)
    feas_sols = (L + np.dot(E, feas_sols_X.T)).T

    return feas_sols_X, feas_sols, Q_I, h, Q, offset


