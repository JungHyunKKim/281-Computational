import numpy as np
import scipy as sp

def firm(calibration, ss, T):

    alpha = calibration['alpha']

    I = sp.sparse.eye(T)
    Ip1 = sp.sparse.diags([np.ones(T-1)], [1], (T, T))
    Im1 = sp.sparse.diags([np.ones(T-1)], [-1], (T, T))
    Z = sp.sparse.csr_matrix((T, T))

    J = {}

    Y = ss['Z'] * ss['K'] ** alpha * ss['L'] ** (1 - alpha)
    r = alpha * Y / ss['K'] - ss['delta']
    w = (1 - alpha) * Y / ss['L']

    # firm block matrices: output
    J['Y'] = {'Z': Y / ss['Z'] * I, 
              'K': alpha * Y / ss['K'] * Im1}
    
    # firm block matrices: real rate
    J['r'] = {'Z': alpha * ss['Y'] / (ss['Z'] * ss['K'])  * I, 
              'K': alpha * (alpha - 1) * ( ss['Y'] / ss['K']**2 ) * Im1}
    
    # firm block matrices: real rate
    J['w'] = {'Z': w / ss['Z'] * I, 
              'K': alpha * w / ss['K'] * Im1}
    
    return J