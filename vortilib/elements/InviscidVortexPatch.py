import numpy as np
    
def omega(X,Y,k=2,bPolarIn=False): 
    if bPolarIn:
        r = X
        theta = Y
    else:
        r = np.sqrt(X ** 2 + Y ** 2)
        theta = np.arctan2(Y,X)
    
    omega_z = (1-r**2)**k
    # r>1
    I = np.abs(r) > 1
    omega_z[I] = 0
    # r=0
    I = np.abs(r) < 1e-14
    omega_z[I] = 1
    return omega_z
