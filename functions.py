from scipy import *
from pylab import *

def max_f(rs,f):
    maxPos = rs[argmax(f)]
    maxVal = max(f)
    return maxPos, maxVal

def energy(rs,m,sigma,f,g,w):
    H = 1.0 - 2.0*m/rs
    Ttot=2.*f**2 - (4.*w**2*f**2)/(H*sigma**2)

    intE = trapz(-rs**2 * sigma * Ttot,x=rs)
    return intE

def charge(rs,m,sigma,f,g,w):
    intQ = trapz(rs**2 * ( 2*f**2/sigma ),x=rs)
    Q = w * intQ
    return Q

def light_rings(rs,m,sigma):
    Veff = ( rs - 2*m ) * sigma**2 / rs**3
    maxPos = argmax(Veff)
    minPos = argmin(Veff)
    if maxPos == 0:
        maxOut = 'No maxima'
    else:
        maxOut = [rs[maxPos],max(Veff)]
    if minPos == len(rs)-1:
        minOut = 'No minima'
    else:
        minOut = [rs[minPos],min(Veff)]
    return minOut,maxOut

def read_parameters():
    pars = genfromtxt('parameters')
    mu = pars[0]
    sigma0 = pars[1]
    w = pars[2]
    a_0 = pars[3]
    r1 = pars[4]
    rmax = pars[5]
    return mu,sigma0,w,a_0,r1,rmax

def plot_functions(rs,m,sigma,f,g):
    plot(rs,f)
    savefig('f.pdf')
    clf()
    plot(rs,g)
    savefig('g.pdf')
