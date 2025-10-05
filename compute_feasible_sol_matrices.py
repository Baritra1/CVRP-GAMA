import numpy as np
from sim_anneal_paths import *
from dwave.system import DWaveSampler

def compute_a(node_num,k,d,q):
    #n^2+n x variables
    #(n-1)log(q) y variables 
    #log(k) slack variables for (13)
    #(n-1)logq slack variables for (15) upper bound
    #Σlogd for (15) lower bound (this is worst case (n+1)logq variables))
    #(n^2-n)(logq+1) slack variables for (14)
    n=node_num+1
    d_copy1=d.copy()
    d_copy1.append(d_copy1[0])
    numbits_d=0
    for i in d_copy1:
        numbits_d+=i.bit_length()
    cols = n*(n-1) + (n-2)*q.bit_length()+k.bit_length()+(n-2)*q.bit_length()+numbits_d+((n-1)*(n-1)-(n-1))*(q.bit_length()+1)
    #11 has n-1 rows
    #corollary to 11 has n rows
    #12 has n-1 rows
    #15 upper bound has n-1 rows
    #15 lower bound has n-1 rows
    #13 has 1 row
    #14 has n^2-n rows
    #made up condition 1 has 1 row
    A = np.zeros(((n-2)+(n-1)+(n-2)+(n-2)+(n-2)+(n-1)*(n-1)-(n-1)+2, cols), dtype=int)
    #A is O(n^2*n^2logq) memory

    #(11)
    for i in range(1,n-1):
        A[i-1][n*i:n*i+n]=1
        A[i][n*i+i]=0 

    #(corollary to 11) x_ii=0
    for i in range(n-1):
        A[(n-2)+i][n*i+i]=1   

    #(12)
    for h in range(1,n-1):
        for i in range(0,n-1):
            if i != h:
                A[(2*n-3)+(h-1)][i*n + h] = 1
        for j in range(1,n):
            if j != h:
                    A[(2*n-3)+(h-1)][h*n + j] = -1

    #(15) upper bound
    for i in range(1,n-1):
        m=1
        for j in range(q.bit_length()):
            A[(3*n-5)+(i-1)][n*(n-1) + (i-1)*q.bit_length()+j]=m
            m*=2
        l=1
        for j in range(q.bit_length()):
            A[(3*n-5)+(i-1)][n*(n-1) + (n-2)*q.bit_length()+k.bit_length()+(i-1)*q.bit_length()+j]=l
            l*=2
    #(15) lower bound
    d_accum=0
    for i in range(1,n-1):
        m=1
        for j in range(q.bit_length()):
            A[(4*n-7)+(i-1)][n*(n-1) + (i-1)*q.bit_length()+j]=m
            m*=2
        l=1
        for j in range(d_copy1[i].bit_length()):
            A[(4*n-7)+(i-1)][n*(n-1) + (n-2)*q.bit_length()+k.bit_length()+(n-2)*q.bit_length()+d_accum+j]=-l
            l*=2
        d_accum+=d_copy1[i].bit_length()
    #14
    accum_14=0
    for i in range(n-1):
        for j in range(n-1):
            if(i==j):
                continue

            l=1
            for h in range(q.bit_length()):
                if(i!=0):
                    A[(5*n-9)+accum_14][n*(n-1)+q.bit_length()*(i-1)+h]=-l
                if(j!=0):
                    A[(5*n-9)+accum_14][n*(n-1)+q.bit_length()*(j-1)+h]=l
                l*=2
            A[(5*n-9)+accum_14][i*n+j]=-(d_copy1[j]+q)
            l=1
            for h in range(q.bit_length()+1):
                A[(5*n-9)+accum_14][n*(n-1) + (n-2)*q.bit_length()+k.bit_length()+(n-2)*q.bit_length()+numbits_d+(accum_14)*(q.bit_length()+1)+h]=-l
                l*=2
            accum_14+=1
    #(13)
    A[(5*n-9)+(n-1)*(n-1)-(n-1)][1:n-1]=1
    l=1
    for i in range(k.bit_length()):
        A[(5*n-9)+(n-1)*(n-1)-(n-1)][n*(n-1)+(n-2)*q.bit_length()+i]=l
        l*=2

    #made up condition 1 - you should not be able to go back to starting spot, only ending spot
    # also, you should not be able to go from start depot to end depot (this is analogous to never leaving)
    for i in range(n-1):
        A[(5*n-9)+(n-1)*(n-1)-(n-1)+1][n*i]=1
        
    A[(5*n-9)+(n-1)*(n-1)-(n-1)+1][n-1]=1

    return A

def compute_b(node_num,k,d,q):
    n=node_num+1
    d_copy2=d.copy()
    d_copy2.append(d_copy2[0])
    b = np.zeros((n-2)+(n-1)+(n-2)+(n-2)+(n-2)+(n-1)*(n-1)-(n-1)+2)
    #(11)
    b[:n-2] = 1
    #(corollary to 11)
    b[n-2:2*n-3]=0
    #(12)
    b[2*n-3:3*n-5]=0
    #(15) upper bound
    b[3*n-5:4*n-7]=q
    #(15) lower bound
    for i in range(1,n-1):
        b[4*n-7+i-1]=d_copy2[i]
    #(14)
    b[5*n-9:(5*n-9)+(n-1)*(n-1)-(n-1)]=-q
    #(13)
    b[(5*n-9)+(n-1)*(n-1)-(n-1)]=k
    #made up condition 1
    b[(5*n-9)+(n-1)*(n-1)-(n-1)+1]=0
    return b

def compute_matrices(k,d,q):
    assert(d[0]==0)
    for i in range(len(d)):
        assert(d[i]<=q)
    assert(k>=len(d))
    A= compute_a(len(d),k,d,q)
    b= compute_b(len(d),k,d,q)
    return A,b

if  __name__ == '__main__':
    A,b = compute_matrices(5,[0,8,12,13,10],24)
    np.set_printoptions(threshold=np.inf)
    # print(A)
    # print(b)
    # arr=[0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0]
    # print(np.dot(A,arr)-b)
    get_feasible(A,b)
