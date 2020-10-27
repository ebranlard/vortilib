"""
Reference:
    [1] E. Branlard - Wind Turbine Aerodynamics and Vorticity Based Method, Springer, 2017
"""
#--- Legacy python 2.7
from __future__ import division
from __future__ import print_function
# --- General
import unittest
import numpy as np

# --------------------------------------------------------------------------------{
def vs_u_raw(CP, Pa, Pb, Gamma, RegFunction=0, RegParam=0):
    """ Induced velocity from a vortex segment on one control point
    See fUi_VortexSegment11_smooth

    CPs : 3 x n

    RegFunction: Regularization function:
                 0: None
                 1: Rankine
                 2: Lamb-Oseen
                 3: Vatistas
                 4: Denominator offset
    """
    DPa = CP.ravel()-Pa.ravel()
    DPb = CP.ravel()-Pb.ravel()
    xa, ya, za = DPa[0], DPa[1], DPa[2]
    xb, yb, zb = DPb[0], DPb[1], DPb[2]

    norm_a      = np.sqrt(xa * xa + ya * ya + za * za)
    norm_b      = np.sqrt(xb * xb + yb * yb + zb * zb)
    denominator = norm_a * norm_b * (norm_a * norm_b + xa * xb + ya * yb + za * zb)

    if (denominator < 1e-17):
        return np.zeros((1,3))
    if (norm_a < 1e-08 or norm_b < 1e-08):
        return np.zeros((1,3))

    crossprod       = np.array([[ya * zb - za * yb, za * xb - xa * zb, xa * yb - ya * xb]])
    # Singular model
    if RegFunction==0:
        Kv = Gamma / (4.0 * np.pi) * (norm_a + norm_b) / denominator
        return Kv * crossprod

    # Regularization models 
    norm2_r0        = (xa - xb)**2 + (ya - yb)**2 + (za - zb)**2
    norm2_crossprod = crossprod[0,0]**2 + crossprod[0,1]**2 + crossprod[0,2]**2
    h2              = norm2_crossprod/norm2_r0 # Orthogonal distance (r1 x r2)/r0
    eps2 = h2/RegParam**2

    if RegFunction==1:
        if (eps2 < 1):
            Kv = eps2
        else:
            Kv = 1.0
    elif RegFunction==2:
        Kv = 1.0 - np.exp(-1.25643 * eps2)
    elif RegFunction==3:
        Kv = eps2 / np.sqrt(1 + eps2**2)
    elif RegFunction==4:
        Kv = 1.0
        denominator = denominator + RegParam**2 * norm2_r0
    Kv = Gamma * Kv / (4.0 * np.pi) * (norm_a + norm_b) / denominator
    return Kv * crossprod

def vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma, RegFunction=0, RegParam=0):
    """ Induced velocity from a vortex segment on several control point
    See fUi_VortexSegment11_smooth

    CPs : n x 3

    RegFunction: Regularization function:
                 0: None
                 1: Rankine
                 2: Lamb-Oseen
                 3: Vatistas
                 4: Denominator offset
    OUTPUTS:
        ux, uy, uz: velocity, shape of Xcp
    """
    Xcp = np.asarray(Xcp)
    shape_in = Xcp.shape
    Xcp = Xcp.ravel()
    Ycp = np.asarray(Ycp).ravel()
    Zcp = np.asarray(Zcp).ravel()
    ux  = np.zeros(Xcp.shape)
    uy  = np.zeros(Xcp.shape)
    uz  = np.zeros(Xcp.shape)
    for i,(x,y,z) in enumerate(zip(Xcp,Ycp,Zcp)):
        CP=np.array([[x,y,z]])
        u = vs_u_raw(CP,Pa,Pb, Gamma, RegFunction, RegParam)
        ux[i] = u[0,0]
        uy[i] = u[0,1]
        uz[i] = u[0,2]
        
    ux = ux.reshape(shape_in)
    uy = uy.reshape(shape_in)
    uz = uz.reshape(shape_in)
    return ux,uy,uz



# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class TestVortexSegment(unittest.TestCase):
    def test_VS_RegFunctions(self):
        import warnings
#         warnings.filterwarnings('error')

        import matplotlib.pyplot as plt
        import matplotlib as mpl
        from pybra.colors import fColrs
        mpl.rcParams.update({'font.size': 15})
        mpl.rcParams['font.size'] = 15

        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.2,4.6)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.98, top=0.96, bottom=0.12, hspace=0.20, wspace=0.20)
        # --- One vortex segment
        for z0 in [1]:
            Pa = np.array([[ 0, 0, -z0]])
            Pb = np.array([[ 0, 0,  z0]])
            # --- test, 0 on singularity
            U  = vs_u_raw(np.array([0,0,0]), Pa, Pb, Gamma = 1, RegFunction = 0, RegParam = 0)
            np.testing.assert_equal(U, np.zeros((1,3)))
            U  = vs_u_raw(Pa, Pa, Pb, Gamma = 1, RegFunction = 0, RegParam = 0)
            np.testing.assert_equal(U, np.zeros((1,3)))
            U  = vs_u_raw(Pb, Pa, Pb, Gamma = 1, RegFunction = 0, RegParam = 0)
            np.testing.assert_equal(U, np.zeros((1,3)))


            # --- Comparison of regularization
            L = np.linalg.norm(Pa-Pb)
            z0=L/2
            print(L)
            for zz0 in [0]:
                Gamma =1
                Epsilon=2
                Epsilon=0.5*L
                Xcp = np.linspace(0,2*L,100)[:]
                Ycp = Xcp*0
                Zcp = Xcp*0+zz0
                U0x, U0y, U0z = vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma=Gamma, RegFunction=0, RegParam=Epsilon)
                U1x, U1y, U1z = vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma=Gamma, RegFunction=1, RegParam=Epsilon)
                U2x, U2y, U2z = vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma=Gamma, RegFunction=2, RegParam=Epsilon)
                U3x, U3y, U3z = vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma=Gamma, RegFunction=3, RegParam=Epsilon)
                U4x, U4y, U4z = vs_u(Xcp, Ycp, Zcp, Pa, Pb, Gamma=Gamma, RegFunction=4, RegParam=Epsilon*1)

                U_th = Gamma*(L/2)/(2*np.pi*Xcp * np.sqrt(Xcp**2 + z0**2))


                ax.plot(Xcp[1:]/L  ,U0y[1:] / Gamma*np.pi*L, '-', color=fColrs(1), label = 'Singular'       )
                ax.plot(Xcp/L  ,U1y / Gamma*np.pi*L, '-.', color=fColrs(2), label = 'Rankine'    )
                ax.plot(Xcp/L  ,U2y / Gamma*np.pi*L, '-', color=fColrs(3), label = 'Lamb-Oseen' )
                ax.plot(Xcp/L  ,U3y / Gamma*np.pi*L, '--', color=fColrs(4), label = 'Vatistas n=2'   )
                ax.plot(Xcp/L  ,U4y / Gamma*np.pi*L, ':', color=fColrs(5), label = 'Denominator')
#                 ax.plot(Xcp/L  ,U_th/ Gamma*np.pi*L, '--',color='k', label = 'Theory'       )
                ax.set_xlabel(r'$\rho/r_0$ [-]')
                ax.set_ylabel(r'$u / (\pi \Gamma  r_0)$ [-]')
                ax.set_xticks(np.arange(0,2.1,0.5))
                ax.set_xticklabels(['0','0.5','1','1.5','2'])
                ax.legend()
                ax.set_ylim([0,1])
                ax.set_xlim([0,2])
                ax.tick_params(direction='in')
        fig.savefig('FilamentRegularization.pdf')
        plt.show()










if __name__ == "__main__":
    unittest.main()


