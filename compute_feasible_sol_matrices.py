import numpy as np
from sim_anneal_paths import *

def compute_a(node_num,k,d,q):
    n=node_num+1
    d_copy1=d.copy()
    d_copy1.append(d_copy1[0])

    #n^2+n x variables
    #log(k) slack variables for (13)
    cols = n*(n-1)+k.bit_length()
    #11 has n-1 rows
    #corollary to 11 has n rows
    #12 has n-1 rows
    #13 has 1 row
    #made up condition 1 has 1 row
    A = np.zeros(((n-2)+(n-1)+(n-2)+2, cols), dtype=int)
    #A is O(n^2*n) memory

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

    #(13)
    A[(3*n-5)][1:n-1]=1
    l=1
    for i in range(k.bit_length()):
        A[(3*n-5)][n*(n-1)+i]=l
        l*=2

    #made up condition 1 - you should not be able to go back to starting spot, only ending spot
    # also, you should not be able to go from start depot to end depot (this is analogous to never leaving)
    for i in range(n-1):
        A[(3*n-5)+1][n*i]=1
        
    A[(3*n-5)+1][n-1]=1

    return A

def compute_b(node_num,k,d,q):
    n=node_num+1
    d_copy2=d.copy()
    d_copy2.append(d_copy2[0])
    b = np.zeros((n-2)+(n-1)+(n-2)+2)
    #(11)
    b[:n-2] = 1
    #(corollary to 11)
    b[n-2:2*n-3]=0
    #(12)
    b[2*n-3:3*n-5]=0

    #(13)
    b[(3*n-5)]=k
    #made up condition 1
    b[(3*n-5)+1]=0
    return b

def compute_matrices(k,d,q):
    assert(d[0]==0)
    sum=0
    for i in range(len(d)):
        assert(d[i]<=q)
    A=compute_a(len(d),k,d,q)
    b=compute_b(len(d),k,d,q)
    return A,b

if  __name__ == '__main__':
    A,b = compute_matrices(57,[0, 7, 2, 4, 3, 5, 7, 6, 7, 2, 0, 3, 0, 1, 0, 3, 3, 7, 5, 7, 5, 7, 1, 7, 6, 3, 8, 3, 0, 3],37)
    np.set_printoptions(threshold=np.inf)
    print(A)
    print(b)
    # arr=[0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0]
    # print(np.dot(A,arr)-b)
    get_feasible(A,b,len([0, 7, 2, 4, 3, 5, 7, 6, 7, 2, 0, 3, 0, 1, 0, 3, 3, 7, 5, 7, 5, 7, 1, 7, 6, 3, 8, 3, 0, 3]),37,[0, 7, 2, 4, 3, 5, 7, 6, 7, 2, 0, 3, 0, 1, 0, 3, 3, 7, 5, 7, 5, 7, 1, 7, 6, 3, 8, 3, 0, 3],10000)
