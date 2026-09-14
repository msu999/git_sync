import numpy as np

K = np.array([
    [1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 2, 2, 3, 4],
    [0, 0, 1, 1, 2, 3, 5],
    [0, 0, 0, 1, 1, 3, 6],
    [0, 0, 0, 0, 1, 2, 5],
    [0, 0, 0, 0, 0, 1, 4],
    [0, 0, 0, 0, 0, 0, 1]
]).transpose()

M = np.array([
    [1, 1, 1, 1, 1, 1],
    [5, 3, 2, 1, 1, 0],
    [10, 4, 1, 2, 0, 0],
    [20, 6, 1, 0, 0, 0],
    [30, 5, 0, 1, 0, 0],
    [60, 1, 0, 0, 0, 0],
    [120, 0, 0, 0, 0, 0]
])

S = np.zeros_like(M)

for (i, r) in enumerate(M):
    res = M[i]

    for j in range(0, i):
        print(f"Iteration {j} in {i} : {K[i, j]*S[j]}")
        res -= K[i, j]*S[j]

    print(f"{i} : res = {res}")

    S[i] = res

print(S)
