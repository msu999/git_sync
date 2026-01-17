import numpy as np

def one(n, i):
    vector = np.zeros(n, dtype=int)
    vector[i] = 1

    return vector

def evalf(f, xs):
    dom_shape = xs.shape

    last_I = dom_shape - np.array([1 for i in dom_shape])
    test_f = f(*xs[*last_I])

    if type(test_f) != np.ndarray:
        res_shape = dom_shape
        res = np.empty(res_shape) 

    else:
        codom_dim = test_f.size
        res_shape = list(dom_shape)
        res = np.empty(res_shape, dtype=np.ndarray)

    #print(f"Test f = {test_f}, Test f size = {codom_dim}")

    print(f"Res shape = {res_shape}")
    

    for (I, x) in np.ndenumerate(xs):
        res[*I] = f(*x)

    return res

def cross(A, B):
    A_shape = A.shape
    B_shape = B.shape
    C_shape = A_shape + B_shape
    print(f"C_shape = {C_shape}")
    C = np.empty(C_shape, dtype=np.ndarray)

    for (I_a, a) in np.ndenumerate(A):
        for (I_b, b) in np.ndenumerate(B):
            if type(a) == np.float64 and type(b) == np.float64:
                C[*I_a, *I_b] = np.array([a, b])

            elif type(b) == np.float64:
                C[*I_a, *I_b] = np.concatenate((a, np.array([b])))
                
            elif type(a) == np.float64:
                C[*I_a, *I_b] = np.concatenate((np.array([a]), b))

            else:
                C[*I_a, *I_b] = np.concatenate((a, b))

    return C

def slicer(A, i, ref_I):
    shape = A.shape  # Shape of array A
    num_vars = len(shape)  # Number of variables
    dim_i = shape[i]  # Dimension of axis i
    test_A_element = A[ref_I]
    result = np.empty(dim_i, dtype=type(test_A_element))  # Array resulting from slicing

    # Slicing
    for j in range(shape[i]):
        I = ref_I + np.zeros(num_vars, dtype=int) 
        I[i] = j
        print(f"I = {I}", end="\n\n")
        result[j] = A[*I]

    return result
