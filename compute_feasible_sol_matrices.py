import numpy as np

def compute_a(node_num):
    #n+1*n+1 variables of x
    #log(n+1) variables for each y, (n+1) y variables
    n=node_num+1
    cols = n*n + n*n.bit_length()
    A = np.zeros((n, cols), dtype=int)

    #(11)
    for i in range(n):
        A[i][n*i:n*i+n]=1
        A[i][n*i+i]=0      
    # x_11-x_12=0
    print(A)

def compute_b(node_num):
    n=node_num+1
    b = np.zeros(n+1)
    #(11)
    b[:n+1] = 1
    #(12)
    b[n+1:n+1+n]=0

if  __name__ == '__main__':
    compute_a(3)
