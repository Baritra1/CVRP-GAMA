import numpy as np

def compute_a(node_num):
    #n+1*n+1 variables of x
    #log(n+1) variables for each y, (n+1) y variables
    n=node_num+1
    A = np.zeros(n*n+n*n.bit_length(),n)

    #(11)
    for i in range(n):
        A[i, i*n:(i+1)*n] = 1
    print(A)

def main():
    compute_a(3)
