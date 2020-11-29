"""
2D Lamb-Oseen vortex

See:
 [1] Chapter 1, p.58, Branlard - Wind turbine aerodynamics and vorticity based methods, Springer 2017


"""

import numpy as np
    
def lo_omega(X,Y,Gamma=1,t=1,nu=1, polarIn=False): 
    """ 
    Vorticity distribution for 2D Lamb-Oseen vortex
    """
    if polarIn:
        r = X
        theta = Y
    else:
        r = np.sqrt(X ** 2 + Y ** 2)
        theta = atan2(Y,X)
    
    omega_z = Gamma/(4*np.pi*nu*t) * (np.exp(- r**2/(4*nu*t)))
    return omega_z
