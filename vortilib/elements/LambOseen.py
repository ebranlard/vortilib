import numpy as np
    
def omega(X,Y,Gamma=1,t=1,nu=1,bPolarIn=False): 
    if bPolarIn:
        r = X
        theta = Y
    else:
        r = np.sqrt(X ** 2 + Y ** 2)
        theta = atan2(Y,X)
    
    omega_z = Gamma/(4*np.pi*nu*t) * (np.exp(- r**2/(4*nu*t)))
    return omega_z
