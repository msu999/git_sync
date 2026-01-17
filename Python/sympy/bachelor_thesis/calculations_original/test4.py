import numpy as np

A = np.array([1, 0, 0])
B = np.array([0, 1, 0])
C = np.array([0, 0, 1])

A_bra = np.diag(A)
B_bra = np.diag(B)
C_bra = np.diag(C)

M = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])

Psi = A + 1j*np.sqrt(2)*B + np.sqrt(2)*C

print(Psi.transpose().conjugate().dot(M))

print(Psi.conjugate().transpose().dot(((A + np.sqrt(2)*B)*B_bra - 1j*C*A_bra)**2).dot(M).dot(Psi))
