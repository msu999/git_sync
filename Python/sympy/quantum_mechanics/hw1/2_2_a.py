from sympy import *

def ip(u, v):
    return u.H.dot(v)

def transfer_matrix(B, C):
    B_1, B_2 = B[0], B[1]
    C_1, C_2 = C[0], C[1]

    return ImmutableMatrix([[ip(B_1, C_1), ip(B_1, C_2)], [ip(B_2, C_1), ip(B_2, C_2)]])

def matrix_in_basis(matrix, transfer_matrix):
    return transfer_matrix.inv() * matrix * transfer_matrix

H = ImmutableMatrix([1, 0])
V = ImmutableMatrix([0, 1])
B_std = (H, V)

D_p = (1/sqrt(2)) * ImmutableMatrix([1, 1])
D_m = (1/sqrt(2)) * ImmutableMatrix([1, -1])
B_D = (D_p, D_m)

CW = (1/sqrt(2)) * ImmutableMatrix([1, -I])
CCW = (1/sqrt(2)) * ImmutableMatrix([1, I])
B_O = (CCW, CW)

T_stdD = transfer_matrix(B_std, B_D)
T_Dstd = T_stdD.inv()
T_stdO = transfer_matrix(B_std, B_O)
T_Ostd = T_stdO.inv()
T_DO = transfer_matrix(B_D, B_O)
T_OD = T_DO.inv()

P_eigen = ImmutableMatrix([[0, 0], [0, 1]])

P_V, P_m, P_CW = P_eigen, P_eigen, P_eigen
P_V_in_D = matrix_in_basis(P_V, T_Dstd)
P_V_in_O = matrix_in_basis(P_V, T_Ostd)
P_m_in_std = matrix_in_basis(P_m, T_stdD)
P_m_in_O = matrix_in_basis(P_m, T_OD)
P_CW_in_std = matrix_in_basis(P_CW, T_stdO)
P_CW_in_D = matrix_in_basis(P_CW, T_DO)

matrices = [[P_V, P_V_in_D, P_V_in_O], [P_m_in_std, P_m, P_m_in_O], [P_CW_in_std, P_CW_in_D, P_CW]]

for r in matrices:
    for matrix in r:
        print(r"\dfrac{1}{2}" + f"{latex(2*matrix.expand().simplify(), mat_delim='(', mat_str='matrix')}\n\n")
