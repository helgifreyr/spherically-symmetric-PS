from scipy import genfromtxt
from pylab import *

dat = genfromtxt('global-data.dat')
w = dat[:,2]
M = dat[:,-1]



# for fake,freq,mass in zip(fake_w,w,M):
#     print fake,freq,mass

plot(w,M,'b.',ms=5)
ylim(0,1.25)
savefig('global-data.pdf')
