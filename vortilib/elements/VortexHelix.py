""" 
Reference:
    [1] E. Branlard - Wind Turbine Aerodynamics and Vorticity Based Method, Springer, 2017, page 477

"""

import numpy as np
import unittest
try:
    from pybra.clean_exceptions import *
except:
    pass
    



def vh_theory_raw_u(r, r0, l, Gamma=1, psih=0, nB=3, bWT=True, bSemi=False, author='wrench'): 
    """ 
    Induced velocity from nB (infinite or semi-infinite) helices, evaluate on a "lifting line"
    Low level function!

    INPUTS:

    NOTE: Equivalent of Matlab fUi_HelixNTheories
    """

    sign = 1
    fact = 1
    if (bWT):
        sign = - 1
    
    if (bSemi):
        fact = 0.5
    
    # Equation 39.12 of [1]
    pexi      = r/r0 * (l + np.sqrt(l**2 + r0**2))/(l + np.sqrt(l**2 + r**2))*np.exp(np.sqrt(l**2 + r**2)/l)/np.exp(np.sqrt(l**2 + r0**2)/l)
    mexi      = 1/pexi
    #Equation 39.14
    C0z       =     ((l**2 + r0**2)/(l**2 + r**2))**(1/4)
    C0r       = 1/l*((l**2 + r0**2)*(l**2 + r**2))**(1/4)
    # Equation 39.17 (Okulov)
    C1r = l/24*((- 2*l**2 - 9*r**2)/(l**2 + r**2)**(3/2) + (2*l**2 + 9*r0**2)/(l**2 + r0**2)**(3/2))
    vr = 0
    vz = 0
    if 'okulov' == (author):
        C1z = C1zWrench
    elif 'wrench' == (author):
        # Equation 39.16
        C1z = l/24*((3*r**2 - 2*l**2)/(l**2 + r**2)**(3/2) + (2*l**2 + 9*r0**2)/(l**2 + r0**2)**(3/2))
    elif 'lerbs' == (author):
        # Equation 39.15
        C1z  = l/2*r0**2/(l**2 + r0**2)**(3/2)
    else:
        raise Exception('No other method')
    
    # Stacking convention is such that upper value corresponds to low r..
    # TODO vectorize me
    if (np.abs(r) < r0):
        try:
            tmp = 1/((mexi*np.exp(-1j*psih))**nB - 1)
        except:
            # Overflow may happen
            #print(mexi,r/r0,nB)
            tmp=0
        vz  =    1/(2*np.pi*l) + 1/(2*np.pi*l)*C0z*np.real( tmp + C1z/nB*np.log(1 + tmp))
        vr  =                  - 1/(2*np.pi*r)*C0r*np.imag( tmp + C1r/nB*np.log(1 + tmp))
    elif (np.abs(r) > r0):
        try:
            tmp = 1/((pexi*np.exp(-1j*psih))**nB - 1)
        except:
            # Overflow may happen (in that case the exponent term is really small)
            #print(mexi,r/r0,nB)
            tmp=0

        vz  =   0              + 1/(2*np.pi*l)*C0z*np.real(-tmp + C1z/nB*np.log(1 + tmp))
        vr  =                  - 1/(2*np.pi*r)*C0r*np.imag( tmp - C1r/nB*np.log(1 + tmp))
    else:
        vr = 0
        vz = 0
    
    vt =      l/r*(1/(2*np.pi*l) - vz)
    vz = sign*vz
    ui = fact*nB*Gamma*np.array([vr,vt,vz])
    return ui


def vh_u(Xcp,Ycp,Zcp,Gamma,R,h,psih=0,nB=3,bWT=True,method='wrench',bSemi=True,polar_out=True):
    """
    Induced velocity by nB azimuthally dstributed helices
    psih : reference azimuthal position of the first helix
    """
    sign = 1
    if (bWT):
        sign = - 1
    l = h/(2*np.pi)

    Xcp=np.asarray(Xcp)
    Ycp=np.asarray(Ycp)
    Zcp=np.asarray(Zcp)


    shape_in=Xcp.shape
    Xcp    = Xcp.ravel()
    Ycp    = Ycp.ravel()
    Zcp    = Zcp.ravel()
    r_cp   = np.sqrt(Xcp**2 + Ycp**2)
    psi_cp = np.arctan2(Ycp,Xcp)
    psih_cp= (-psih + psi_cp) - sign * Zcp / l

    ur  = np.zeros(r_cp.shape)
    ut  = np.zeros(r_cp.shape)
    uz  = np.zeros(r_cp.shape)
    # ---- Loop on all control points to find velocity
    for i,(rc,psihc,zc) in enumerate(zip(r_cp,psih_cp,Zcp)):
        #      ui[i,:-1] = ui(i,:) + fUi_HelixNTheories(Algo.Helix.Method(np.arange(1,end() - 1+1)),Gamma,rc,R,l,,nB,bWT,bHalf)
        u=vh_theory_raw_u(r=rc, r0=R, l=l, Gamma=Gamma, psih=psihc, nB=nB, bWT=bWT , bSemi=bSemi, author=method)
        ur[i]+=u[0]
        ut[i]+=u[1]
        uz[i]+=u[2]

    ur=ur.reshape(shape_in)
    ut=ut.reshape(shape_in)
    uz=uz.reshape(shape_in)
    if polar_out:
        return ur,ut,uz
    else:
        raise NotImplementedError()


class TestHelix(unittest.TestCase):

    def Params(self,a=1/3,U0=10,Lambda=6,R=100,nB=3):
        h         = 2*np.pi*R*(1-a)/Lambda
        CT        = 4*a*(1-a)
        Gamma     = CT*np.pi*R*U0/(nB*Lambda) # NOTE assume aprime = 0, large tip-sped ratio 
        return  a,U0,R,Lambda,nB,h,CT,Gamma


    def test_VH_LiftingLine(self):
        # TODO radial induction is zero on the lifting line due to "inherent" infinite assumption from theory
        import warnings
        warnings.filterwarnings('error')

        a,U0,R,Lambda,nB,h,CT,Gamma_B=self.Params(nB=10)
        method = 'wrench'
        bWT   = True # Convention WT or Propeller
        bSemi = True # infinite helix or at the rotor

        # ---- On the lifting line, close to the root
        # The axial induction is close to the one for an infinite number of blades 
        # The larger B the more it is constant along the span (ie closer to a vortex cylinder)
        for nB in [1,3,10]:
            Gamma_B   = CT*np.pi*R*U0/(nB*Lambda) # NOTE assume aprime = 0, large tip-sped ratio
            gamma_t = - nB*Gamma_B/h # Vortex cylinder (if infinite number of blades)
            ur,ut,uz = vh_u([0.01],[0],[0],Gamma_B,R,h,psih=0,nB=nB,bWT=bWT,method=method,bSemi=bSemi)
            np.testing.assert_almost_equal(uz,[-U0*a], 5)
            np.testing.assert_almost_equal(uz,[gamma_t/2],5)
            if nB>9:
                ur,ut,uz = vh_u([0.5*R],[0],[0],Gamma_B,R,h,psih=0,nB=nB,bWT=bWT,method=method,bSemi=bSemi)
                np.testing.assert_almost_equal(uz,[-U0*a], 5)

        # ---- On the lifting line outside
        # The axial induction should be zero when r>1.5R (approx)
        for nB in [3,10]:
            Gamma_B   = CT*np.pi*R*U0/(nB*Lambda) # NOTE assume aprime = 0, large tip-sped ratio
            gamma_t = - nB*Gamma_B/h # Vortex cylinder (if infinite number of blades)
            ur,ut,uz = vh_u([R*1.5],[0],[0],Gamma_B,R,h,psih=0,nB=nB,bWT=bWT,method=method,bSemi=bSemi)
            np.testing.assert_almost_equal(uz,[0], 4)


        # --- Velocity on lifting line should be half the one of infinite helix
        vr      = np.arange(0.1*R,2*R,0.2*R)
        nB      = 3
        Gamma_B = CT*np.pi*R*U0/(nB*Lambda)  # NOTE assume aprime = 0, large tip-sped ratio
        ur ,ut ,uz  = vh_u([vr],[vr*0],[vr*0],Gamma_B,R,h,psih=0,nB=nB,bWT=bWT,method=method,bSemi=True)
        ur2,ut2,uz2 = vh_u([vr],[vr*0],[vr*0],Gamma_B,R,h,psih=0,nB=nB,bWT=bWT,method=method,bSemi=False)
        np.testing.assert_almost_equal(2*uz,uz2, 4)
        np.testing.assert_almost_equal(2*ut,ut2, 4)



if __name__ == '__main__':
    unittest.main()

    
# def fUi_HelixN(nB = None,v1 = None,v2 = None,v3 = None,vGamma = None,vR = None,vh = None,vpsi = None,bGridIn = None,bPolarIn = None,bPolarOut = None,bComputeGrad = None,Algo = None): 
#     # N is a nasty argument
# # h is my ie h=2pi r tan(eps)
#     
#     # TODO fUi_HelixN(1,0,0,-1,4,1,0.1,0,0,0,0,0,Algo) > returns NaN
#     
#     bWT = Algo.Helix.bWT
#     bHalf = not Algo.Helix.bInf 
#     sign = 1
#     if (bWT):
#         sign = - 1
#     
#     nHelix = len(vh)
#     if (not (Algo.Helix.Method(end()) == 'N') ):
#         print('Computing %d Helix using for loop' % (nHelix))
#         tic()
#         ui = 0
#         for i in np.arange(1,nHelix+1).reshape(-1):
#             uiw,grad,Xcp,Ycp,Zcp,Helix[i-1] = fUi_Helix1(v1,v2,v3,vGamma(i),vR(i),vh(i),vpsi(i),bGridIn,bPolarIn,bPolarOut,bComputeGrad,Algo)
#             ui = uiw + ui
#             if (np.mod(i,3) == 1):
#                 print('.' % ())
#         toc()
#     else:
#         #     if(length(vh)==N && N~=1)
# #         error('I guess I should implement a c code tht does all helix at once');
# #     else
#         print('Computing %d systems of %d helix using theory:' % (nHelix,nB))
#         tic()
#         if (nargout > 4):
#             theta = np.linspace(0,Algo.nRev * 2 * pi,Algo.nRev * Algo.nPhi)
#             theta = theta * sign
#             #we have to compute the Helixes
#             for i in np.arange(1,nHelix+1).reshape(-1):
#                 vpsi_helix = np.arange(0,(nB - 1) * 2 * pi / nB + vpsi(i)+1)
#                 l = vh(i) / (2 * pi)
#                 R = vR(i)
#                 #helix coordinates
#                 for b in np.arange(1,nB+1).reshape(-1):
#                     xh = R * np.cos(theta + vpsi_helix(b))
#                     yh = R * np.sin(theta + vpsi_helix(b))
#                     zh = np.abs(l * theta)
#                     Helix[i-1][b-1] = np.array([np.transpose(xh),np.transpose(yh),np.transpose(zh)])
#         Xcp_flat,Ycp_flat,Zcp_flat,ncp,Xcp,Ycp,Zcp = fgetControlPoints(v1,v2,v3,bPolarIn,bGridIn)
#         if bComputeGrad:
#             # TODO
#             grad = np.zeros((ncp,9))
#         else:
#             grad = []
#         ui = np.zeros((ncp,3))
#         if 'okulovbN' == (Algo.Helix.Method):
#             for i in np.arange(1,nHelix+1).reshape(-1):
#                 l = vh(i) / (2 * pi)
#                 psih = vpsi(i)
#                 R = vR(i)
#                 Gamma = vGamma(i)
#                 for i in np.arange(1,ncp+1).reshape(-1):
#                     zc = Zcp_flat(i)
#                     rc = np.sqrt(Ycp_flat(i) ** 2 + Xcp_flat(i) ** 2)
#                     psic = atan2(Ycp_flat(i),Xcp_flat(i))
#                     ui[i,:-1] = ui(i,:) + Gamma * fUiOkulovN(rc,R,l,(- psih + psic) - sign * zc / l,nB,bWT,bHalf)
#         else:
#             for i in np.arange(1,nHelix+1).reshape(-1):
#                 l = vh(i) / (2 * pi)
#                 psih = vpsi(i)
#                 R = vR(i)
#                 Gamma = vGamma(i)
#                 for i in np.arange(1,ncp+1).reshape(-1):
#                     zc = Zcp_flat(i)
#                     rc = np.sqrt(Ycp_flat(i) ** 2 + Xcp_flat(i) ** 2)
#                     psic = atan2(Ycp_flat(i),Xcp_flat(i))
#                     ui[i,:-1] = ui(i,:) + fUi_HelixNTheories(Algo.Helix.Method(np.arange(1,end() - 1+1)),Gamma,rc,R,l,(- psih + psic) - sign * zc / l,nB,bWT,bHalf)
#         ## TODO understand this
#         in_ = np.isnan(ui(:,1))
#         ui[in_,1-1] = 0
#         in_ = np.isnan(ui(:,2))
#         ui[in_,2-1] = 0
#         if (not bPolarOut ):
#             phi = atan2(Ycp_flat,Xcp_flat)
#             urad = ui(:,1)
#             upsi = ui(:,2)
#             ui[:,1-1] = np.multiply(upsi,np.sin(phi)) + np.multiply(urad,np.cos(phi))
#             ui[:,2-1] = np.multiply(upsi,np.cos(phi)) - np.multiply(urad,np.sin(phi))
#         #  At the end we reshape
#         if (bGridIn):
#             ui = fReshape(ui,len(v1),len(v2),len(v3),3)
#         toc()
#     
#     return ui,grad,Xcp,Ycp,Zcp,Helix
#     
#     ##########################################################################
# ###
# ##########################################################################
#     
# def fUiOkulovN(r = None,a = None,l = None,psih = None,nB = None,bWT = None,bHalf = None): 
#     sign = 1
#     fact = 1
#     if (bWT):
#         sign = - 1
#     
#     if (bHalf):
#         fact = 0.5
#     
#     Ar = ((l ** 2 + r ** 2) * (l ** 2 + a ** 2)) ** (1 / 4) / l
#     A = (l ** 2 + a ** 2) ** (1 / 4) / (l ** 2 + r ** 2) ** (1 / 4)
#     C = l / 24 * ((- 2 * l ** 2 - 9 * r ** 2) / (l ** 2 + r ** 2) ** (3 / 2) + (2 * l ** 2 + 9 * a ** 2) / (l ** 2 + a ** 2) ** (3 / 2))
#     B = l / 24 * ((3 * r ** 2 - 2 * l ** 2) / (l ** 2 + r ** 2) ** (3 / 2) + (2 * l ** 2 + 9 * a ** 2) / (l ** 2 + a ** 2) ** (3 / 2))
#     vr = 0
#     vz = 0
#     for iB in np.arange(1,nB+1).reshape(-1):
#         t = psih + (iB - 1) * 2 * pi / nB
#         if (np.abs(r) < a):
#             vr = vr + - 1 / (2 * pi * r) * Ar * imag(fSm(r / l,a / l,t) + C * fSl(r / l,a / l,t))
#             vz = vz + 1 + A * real(fSm(r / l,a / l,t) + B * fSl(r / l,a / l,t))
#         else:
#             vr = vr + - 1 / (2 * pi * r) * Ar * imag(fSm(a / l,r / l,t) - C * fSl(a / l,r / l,t))
#             vz = vz + 0 + A * real(- fSm(a / l,r / l,t) + B * fSl(a / l,r / l,t))
#     
#     vz = fact * vz * 1 / (2 * pi * l)
#     vt = l / r * (nB * fact * 1 / (2 * pi * l) - vz)
#     vz = sign * vz
#     ui = np.array([vr,vt,vz])
#     return ui
#     
#     ##########################################################################
# ###
# ##########################################################################
#     
# def fSl(x = None,y = None,t = None): 
#     Sl = - np.log(1 - np.exp((np.log(x / y * (np.sqrt(1 + y ** 2) + 1) / (np.sqrt(1 + x ** 2) + 1)) + np.sqrt(1 + x ** 2) - np.sqrt(1 + y ** 2)) + 1j * t))
#     return Sl
#     
#     # fSl(a/l,r/l,t)
# ##########################################################################
# ###
# ##########################################################################
#     
# def fSm(x = None,y = None,t = None): 
#     # Sm=exp(1i*t)/( exp( -(  log(x/y * (sqrt(1+y^2) +1 )/(sqrt(1+x^2)+1)  ) +sqrt(1+x^2)-sqrt(1+y^2) ) )  - exp(1i*t)  );
#     Sm = 1 / (np.exp(- (np.log(x / y * (np.sqrt(1 + y ** 2) + 1) / (np.sqrt(1 + x ** 2) + 1)) + np.sqrt(1 + x ** 2) - np.sqrt(1 + y ** 2))) * np.exp(- 1j * t) - 1)
#     return Sm
#     
#     return ui,grad,Xcp,Ycp,Zcp,Helix

    
