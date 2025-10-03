import numpy as np
from sim_anneal_paths import *
from dwave.system import DWaveSampler

def compute_a(node_num,k,d,q):

    assert(d[0]==0)
    assert(node_num==len(d))
    #n+1*n+1 variables of x
    #log(q) variables for each y, (n+1) y variables
    #log(k) slack variables for (13)
    #(n+1)logq slack variables for (15) upper bound
    #numbits_d slack variables for (15) lower bound, worst case (n+1)logq (chances are youll need <nlogq, reason behind dynamic assignment of variables)
    #n(n-1)(log(q)+1)/2 slack variables for (14)
    # i have a suspicion that i have off by ones here
    n=node_num+1
    d.append(d[0])
    numbits_d=0
    for i in d:
        numbits_d+=i.bit_length()
    # we might be able to save some variables in secon n*q.bit_length()
    cols = n*(n-1) + n*q.bit_length()+k.bit_length()+n*q.bit_length()+numbits_d+(n*(n-1))*(q.bit_length()+1)
    #11 has n-1 rows
    #corollary to 11 has n rows
    #12 has n-1 rows
    #15 upper bound has n rows
    #15 lower bound has n rows
    #13 has 1 row
    #14 has n(n-1) rows
    #made up condition 1 has 1 row
    A = np.zeros(((n-2)+(n-2)+n+n+n-1+(n-1)*(n-1), cols), dtype=int)

    #(11)
    for i in range(1,n-1):
        A[i-1][n*i:n*i+n]=1
        A[i][n*i+i]=0 

    #(corollary to 11) x_ii=0
    for i in range(n-1):
        A[n-2+i][n*i+i]=1   

    #(12)
    for h in range(1,n-1):
        for i in range(0,n-1):
            if i != h:
                A[2*n-4+h][i*n + h] = 1
        for j in range(1,n):
            if j != h:
                    A[2*n-4+h][h*n + j] = -1

    #(15) upper bound
    for i in range(n):
        m=1
        for j in range(q.bit_length()):
            A[3*n-5+i][n*(n-1) + i*q.bit_length()+j]=m
            m*=2
        l=1
        for j in range(q.bit_length()):
            A[3*n-5+i][n*(n-1) + n*q.bit_length()+k.bit_length()+i*q.bit_length()+j]=l
            l*=2
    #(15) lower bound
    d_accum=0
    for i in range(n):
        m=1
        if d[i]==0:
            continue
        for j in range(q.bit_length()):
            A[4*n-6+i][n*(n-1) + i*q.bit_length()+j]=m
            m*=2
        l=1
        for j in range(d[i].bit_length()):
            A[4*n-6+i][n*(n-1) + n*q.bit_length()+k.bit_length()+n*q.bit_length()+d_accum+j]=-l
            l*=2
        d_accum+=d[i].bit_length()
    #14
    accum_14=0
    for i in range(n-1):
        for j in range(n):
            if(i==j):
                continue
            l=1
            for h in range(q.bit_length()):
                A[5*n-7+accum_14][n*(n-1)+q.bit_length()*i+h]=-l
                A[5*n-7+accum_14][n*(n-1)+q.bit_length()*j+h]=l
                l*=2
            A[5*n-7+accum_14][i*n+j]=-(d[j]+q)
            l=1
            for h in range(q.bit_length()+1):
                A[5*n-7+accum_14][n*(n-1) + n*q.bit_length()+k.bit_length()+n*q.bit_length()+numbits_d+(accum_14)*(q.bit_length()+1)+h]=-l
                l*=2
            accum_14+=1
            
    #(13)
    A[5*n-7+(n-1)*(n-1)][1:n-1]=1
    l=1
    for i in range(k.bit_length()):
        A[5*n-7+(n-1)*(n-1)][n*(n-1)+n*q.bit_length()+i]=l
        l*=2

    #made up condition 1 - you should not be able to go back to starting spot, only ending spot
    # also, you should not be able to go from start depot to end depot (this is analogous to never leaving)
    for i in range(n-1):
        A[5*n-6+(n-1)*(n-1)][n*i]=1
        
    A[5*n-6+(n-1)*(n-1)][n-1]=1

    return A

def compute_b(node_num,k,d,q):
    
    assert(d[0]==0)
    assert(node_num==len(d))

    n=node_num+1
    d.append(d[0])
    b = np.zeros(5*n-4+n*(n-2))
    #(11)
    b[:n-2] = 1
    #(corollary to 11)
    b[n-2:2*n-3]=0
    #(12)
    b[2*n-3:3*n-5]=0
    #(15) upper bound
    b[3*n-5:4*n-5]=q
    #(15) lower bound
    for i in range(1,n):
        b[4*n-5+i-1]=d[i]
    #(14)
    b[5*n-7:5*n+n*(n-2)-6]=-q
    #(13)
    b[5*n-6+n*(n-2)]=k
    #made up condition 1
    b[5*n-6+n*(n-2)+1]=0
    return b

if  __name__ == '__main__':
    A = compute_a(3,5,[0,20,19],24)
    b = compute_b(3,5,[0,20,19],24)
    np.set_printoptions(threshold=np.inf)
    print(A)
    print(b)
    get_feasible(A,b)
