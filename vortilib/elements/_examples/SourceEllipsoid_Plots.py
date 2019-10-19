from vortilib.elements.SourceEllipsoid import *
import numpy as np
import matplotlib.pyplot as plt
try:
    from pybra.curves import streamQuiver
    from pybra.tictoc import Timer
    from pybra.figure import * 
    setFigureFont(15)
    setFigurePath('./')
except:
    def streamQuiver(*args,**kwargs):
        pass

# --- Parameters
minSpeed=0
maxSpeed=1.20


nx=200
ny=nx+1
# Ellipse parameters
a  = 1
b  = 0.5*a
U0 = 10


# --- Velocity field on grid
vx = np.linspace(-2*a,2*a,nx)
vy = np.linspace(-4*b   ,4*b,ny)
X,Y = np.meshgrid(vx, vy)
# with Timer('Numerical'):
U1,V1 = ser_u_numerical(X,Y,vx,vy,U0,a,b)
# with Timer('Analytical'):
U,V   = ser_u(X,Y,U0,a,b)

# PSI = ser_psi_elliptic(MU,ZETA,U0,a,e)
# PSI = ser_psi(X,Y,X*0,U0,a,b)
# PHI = ser_phi_elliptic(MU,ZETA,U0,a,e)


# --- Plot
Utot     = U+U0
Speed = np.sqrt((Utot**2+V**2))/U0
bInEllipse=(X**2/a**2+Y**2/b**2)<1
Speed[bInEllipse]=np.nan

xe,ye=ellipse_coord(a,b)
fig,ax = plt.subplots(1,1)
# im = ax.pcolormesh(X, Y, Speed, vmin=minSpeed, vmax=maxSpeed)
im = ax.contourf(X, Y, Speed, levels=np.linspace(minSpeed,maxSpeed,25), vmin=minSpeed, vmax=maxSpeed)
cb=fig.colorbar(im)
yseed=np.linspace(np.min(vy)*0.85,np.max(vy)*0.85,8)
start=np.array([yseed*0-2*a*0.9,yseed])
sp=ax.streamplot(vx,vy,Utot,V,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
qv=streamQuiver(ax,sp,n=7,scale=40,angles='xy')
ax.plot(xe,ye,'k',lw=3)
ax.set_ylim([-4*b,4*b])
ax.set_xlim([-2*a,2*a])
ax.set_xlabel('x/a [-]')
ax.set_ylabel('r/a [-]')
ax.set_aspect('equal','box')
ax.set_title('Source Ellipsoid Streamlines')


try:
    export2pdf()
except:
    pass

plt.show()




