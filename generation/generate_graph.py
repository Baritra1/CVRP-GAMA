import numpy as np


#
def generate_graph(n,clamp=1):
    """
    :param clamp: Constrains maximum value of d[i], use to greatly increase # of feasible solutions
    """
    length = n

    q = np.random.randint(1, 101)

    d = np.zeros(length, dtype=int)
    if length > 1:
        d[1:] = np.random.randint(0, (q + 1)/clamp, size=length - 1)
    d=d.tolist()

    k = np.random.randint(length, length + 50)  # k >= len(d) 
    print("d:", d)
    print("q:", q)
    print("k:", k)
    return q,d,k
