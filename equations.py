from scipy import *
from scipy.integrate import odeint
from scipy.optimize import fsolve
from scipy.integrate import trapz

from functions import *

def function(U,r,mu,sigma0,w):
    u1,u2,u3,u4 = U
    # u1 = m, u2 = sigma, u3 = f, u4 = g
    H = 1.0-2.0*u1/r
    du2dr = mu*mu*r*(u4*u4 + u3*u3/(H*H * u2*u2))*u2
    du3dr = (u4*(-(H * mu*mu * u2*u2) + w*w))/w
    du1dr = r*r * ( mu*mu * (u4*u4 * H + u3*u3 / ( H * u2*u2 ) ) \
            + ( du3dr - u4*w )**2 / (u2*u2) )/2.0
    dH = -2.0*du1dr/r + 2.0*u1/(r*r)
    du4dr = u4 * ( -(dH/H) - 2.0/r - du2dr/u2 ) - u3*w / ( H*H * u2*u2 ) 
    # print '%12.11f'%r,'%12.11f'%u1,'%12.11f'%u2,'%12.11f'%u3,'%12.11f'%u4
    return [du1dr, du2dr, du3dr, du4dr]


def odeShooting(a,pars,r2,printStuff=0):
    mu,sigma0,w,a_0,r1,rmax = pars
    r2s = concatenate((array([5.0,10.0]),arange(11,rmax+1)))
    # u1 = m, u2 = sigma, u3 = f, u4 = g
    a = a[0]

    m3 = a*a * mu*mu / ( 6.0 * sigma0*sigma0 )
    sigma2 = a*a * mu*mu  / ( 2.0 * sigma0*sigma0 )
    f2 = a*( mu*mu - w*w/( sigma0*sigma0) )/6.0
    g1 = -a*w/( 3.0 * sigma0*sigma0 )

    # u1 = m, u2 = sigma, u3 = f, u4 = g
    u1_0 = m3 * r1**3
    u2_0 = sigma0 + sigma2 * r1*r1
    u3_0 = a + f2 * r1*r1
    u3_1 = 0.0
    u4_0 = -g1 * r1
    # print u1_0, u2_0, u3_0, u4_0

    rs = arange(r1,r2,0.001)

    U = \
    odeint(function,[u1_0,u2_0,u3_0,u4_0],rs,args=(mu,sigma0,w),atol=1e-12,rtol=1e-12,h0=1e-3)
    # print U
    u1 = U[:,0]
    u2 = U[:,1]
    u3 = U[:,2]
    u4 = U[:,3]
    if printStuff==0:
        return u3[-1]-u3_1
    else:
        return rs,u1,u2,u3,u4
