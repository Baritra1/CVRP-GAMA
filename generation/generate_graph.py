import numpy as np

length = 11

q = np.random.randint(1, 101)

d = np.zeros(length, dtype=int)
if length > 1:
    d[1:] = np.random.randint(0, (q + 1)/4, size=length - 1)

k = np.random.randint(length, length + 50)  # k >= len(d)

print("d:", d)
print("q:", q)
print("k:", k)