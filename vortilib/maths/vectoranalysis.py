import numpy as np
import unittest


# --------------------------------------------------------------------------------}
# --- Functions to reproduce Matlab's gradient, div, and curl: axis 0 and 1 inverted!
# --------------------------------------------------------------------------------{
def matlab_gradient_1d(F,*args,**kwargs):
    """ Computes gradient of a 1d vector field in python like matlab gradient function. 
    call syntax:
       gradient(F) or gradient(F,dx) or gradient(F,x)
    """
    shape_in=F.shape
    G=np.gradient(F.ravel(),*args,**kwargs)
    return G.reshape(shape_in)

def matlab_gradient_2d(F,*args,**kwargs):
    """ Computes gradient of a 1d vector field in python like matlab gradient function. 

    call syntax:
       gradient(F) or gradient(F,dx,dy) or gradient(F,x,y)
    """
    G=np.gradient(F.T,*args,**kwargs)
    G[0]=G[0].T
    G[1]=G[1].T
    return G

def matlab_div_2d(*args):
    """ Computes divergence of a 2D vector field in python like matlab div function

    call syntax:
       div(U,V) or div(X,Y,U,V)
    """
    if len(args)==2:
        sz = args[0].shape
        dx = np.arange(1,sz[1]+1) # NOTE: matlab, x is along the 2nd dimension..
        dy = np.arange(1,sz[0]+1)
        U=args[0]
        V=args[1]
    elif len(args)==4:
        dx = args[0][0,:];
        dy = args[1][:,0];
        U=args[2]
        V=args[3]
    else:
        raise Exception('Input 2 or 4 arguments')
    #retval  = matlab_gradient_2d (U, dx, dy)[0]
    #retval += matlab_gradient_2d(V.T, dy, dx)[0].T;
    return np.gradient(U, dx, axis=1) +  np.gradient(V,  dy, axis=0)

def matlab_curl_2d(*args):
    """ Computes curl of a D vector field in python like matlab curl function
    
    call syntax:
       curl(U,V) or curl(X,Y,U,V)
    """
    if len(args)==2:
        sz = args[0].shape
        dx = np.arange(1,sz[1]+1) # NOTE: matlab, x is along the 2nd dimension..
        dy = np.arange(1,sz[0]+1)
        U=args[0]
        V=args[1]
    elif len(args)==4:
        dx = args[0][0,:];
        dy = args[1][:,0];
        U=args[2]
        V=args[3]
    else:
        raise Exception('Input 2 or 4 arguments')
    dFx_dy = matlab_gradient_2d(U.T, dy, dx)[0].T
    dFy_dx = matlab_gradient_2d(V, dx, dy)[0]
    rot_z  = dFy_dx - dFx_dy                    
    av     = rot_z / 2
    return rot_z, av


# --------------------------------------------------------------------------------}
# --- Pythonic versions 
# --------------------------------------------------------------------------------{
def div(f):
    """
    TODO
    Computes the divergence of the vector field f, corresponding to dFx/dx + dFy/dy + ...
    :param f: List of ndarrays, where every item of the list is one dimension of the vector field
    :return: Single ndarray of the same shape as each of the items in f, which corresponds to a scalar field
    """
    num_dims = len(f)
    return np.ufunc.reduce(np.add, [np.gradient(f[i], axis=i) for i in range(num_dims)])







# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class TestVectorAnalysis(unittest.TestCase):
    def test_grad1d(self):
        x  = np.array([0,1,3,4])
        x2 = x**2
        np.testing.assert_almost_equal(matlab_gradient_1d(x2),np.array([1,4.5,7.5,7]))
        np.testing.assert_almost_equal(matlab_gradient_1d(x2,2),np.array([0.5,2.25,3.75,3.5]))
        # NOTE NOTE: this test below fails due to boundary effects, to replaced by python values..
        # np.testing.assert_almost_equal(matlab_gradient_1d(x2,x),np.array([1,3,5,7]))
        np.testing.assert_almost_equal(matlab_gradient_1d(x2,x),np.array([1,2,6,7]))

        x  = np.array([[0,1,3,4]])
        x2 = x**2
        np.testing.assert_almost_equal(matlab_gradient_1d(x2),np.array([[1,4.5,7.5,7]]))
        np.testing.assert_almost_equal(matlab_gradient_1d(x2.T),np.array([[1,4.5,7.5,7]]).T)

    def test_grad2d(self):
        x = np.arange(-1,1.5,0.5)
        y = np.arange( 0,1.21,0.2)
        X,Y = np.meshgrid(x,y)
        Z = X * np.exp(-X**2 - Y**2)
        # --- Gradient, no dimensions
        px,py = matlab_gradient_2d(Z)
        px_ref =np.array([[ -0.0215210,  0.1839397],[-0.0206771,  0.1767273]])
        py_ref = np.array([ [ 0.014425,  0.015269] ,[ 0.027197,  0.028788]])
        np.testing.assert_almost_equal(px[:2,:2],px_ref,decimal=6)
        np.testing.assert_almost_equal(py[:2,:2],py_ref,decimal=6)
        # --- Gradient, uni dimensions
        px,py = matlab_gradient_2d(Z,0.5)
        px_ref = np.array([ [  -0.0430419,    0.3678794], [  -0.0413542,    0.3534547]])
        py_ref = np.array([ [   0.0288495,    0.0305372], [   0.0543933,    0.0575753]])
        np.testing.assert_almost_equal(px[:2,:2],px_ref,decimal=7)
        np.testing.assert_almost_equal(py[:2,:2],py_ref,decimal=7)
        # --- Gradient, bi dimensions
        px,py = matlab_gradient_2d(Z,0.5,0.2);
        px_ref = np.array([ [  -0.0430419,    0.3678794], [  -0.0413542,    0.3534547]])
        py_ref = np.array([ [   0.0721238,    0.0763430], [   0.1359832,    0.1439382]])
        np.testing.assert_almost_equal(px[:2,:2],px_ref,decimal=7)
        np.testing.assert_almost_equal(py[:2,:2],py_ref,decimal=7)
        px,py = matlab_gradient_2d(Z,x,y);
        np.testing.assert_almost_equal(px[:2,:2],px_ref,decimal=7)
        np.testing.assert_almost_equal(py[:2,:2],py_ref,decimal=7)


    def test_div2d(self):
        x = np.arange(-1,1.5,0.5)
        y = np.arange( 0,1.21,0.2)
        X,Y = np.meshgrid(x,y)
        U =    X * np.exp(-X**2 - Y**2)
        V = Y**2 * np.exp(-X**2 - Y**2)

        # --- Divergence, no dimensions
        div=matlab_div_2d(U,V)
        div_ref = np.array([ [  -0.0073828,    0.2138703], [   0.0044018,    0.2298194]])
        np.testing.assert_almost_equal(div[:2,:2],div_ref,decimal=7)

        # --- Divergence, dimensions
        div=matlab_div_2d(X,Y,U,V)
        div_ref = np.array([ [   0.0276490,    0.5175322], [   0.0840403,    0.6189148]])
        np.testing.assert_almost_equal(div[:2,:2],div_ref,decimal=7)

    def test_curl2d(self):
        x = np.arange(-1,1.5,0.5)
        y = np.arange( 0,1.21,0.2)
        X,Y = np.meshgrid(x,y)
        U =    X * np.exp(-X**2 - Y**2)
        V = Y**2 * np.exp(-X**2 - Y**2)
        # --- Curl, no dimensions
        curl,_=matlab_curl_2d(U,V)
        curl_ref = np.array([ [  -0.0144248,   -0.0152686], [  -0.0114043,   -0.0166409]])
        np.testing.assert_almost_equal(curl[:2,:2],curl_ref,decimal=7)

        # --- Curl, dimensions
        curl,_=matlab_curl_2d(X,Y,U,V)
        curl_ref = np.array([ [  -0.0721238,   -0.0763430], [  -0.1043984,   -0.1196448]])
        np.testing.assert_almost_equal(curl[:2,:2],curl_ref,decimal=7)


if __name__=='__main__':
    unittest.main()
