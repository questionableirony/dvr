import numpy as np
import math, cmath
import matplotlib.pyplot as plt
import time
from scipy import linalg as sclinalg
#from mpl_toolkits.mplot3d import Axes3D
#import scipy as sp #.signal#.lfilter
#from scipy import signal
#from numpy import matrix
#from numpy import linalg

N = 100

X = np.zeros((N,N))

# Hermite polynomial
X[0,1] = np.sqrt(0.5)
X[N-1,N-2] = np.sqrt(0.5*float(N-1))
for i in range(1,N-1):
    X[i,i+1] = np.sqrt(0.5*float(i+1))
    X[i,i-1] = np.sqrt(0.5*float(i))

x, U = np.linalg.eig(X)
Ucc = np.transpose(np.conjugate(U))
# note the x ostensibly covers -inf to inf
# check: np.linalg.inv(U) - np.transpose(np.conjugate(U))   # should be zero
#print np.dot(np.transpose(np.conjugate(U)),np.dot(X,U))    
#print x                                                    # should be the same (along diagonal)!

hbar = 1
mass = 1
omega = 1
VDVR = np.diag(0.5*mass*omega**2*x**2)

Tphi = np.zeros((N,N))
Tphi[0,0] = -0.5*hbar**2/mass * -0.5
Tphi[0,2] = -0.5*hbar**2/mass * 0.5*np.sqrt(2.0)
Tphi[1,1] = -0.5*hbar**2/mass * -1.5
Tphi[1,3] = -0.5*hbar**2/mass * 0.5*np.sqrt(3.0*2.0)
Tphi[N-2,N-2] = -0.5*hbar**2/mass * -(float(N-2)+0.5)
Tphi[N-2,N-4] = -0.5*hbar**2/mass * 0.5*np.sqrt(float(N-2)*float(N-3))
Tphi[N-1,N-1] = -0.5*hbar**2/mass * -(float(N-1)+0.5)
Tphi[N-1,N-3] = -0.5*hbar**2/mass * 0.5*np.sqrt(float(N-1)*float(N-3))
for i in range(2,N-2):
    Tphi[i,i] = -0.5*hbar**2/mass * -(float(i)+0.5)
    Tphi[i,i+2] = -0.5*hbar**2/mass * 0.5*np.sqrt(float(i+2)*float(i+1))
    Tphi[i,i-2] = -0.5*hbar**2/mass * 0.5*np.sqrt(float(i)*float(i-1))
# Tphi appears to be tridiagonalish

TDVR = np.dot(Ucc,np.dot(Tphi,U))

HDVR = TDVR + VDVR

E, eigstates = np.linalg.eig(HDVR)

print np.sort(E)

print hbar*omega*(np.array(range(N))+0.5)

raw_input("")

alpha = 1.0
xbeta = 0.0
pbeta = 0.0
Phi = (np.pi/2.0/alpha)**0.25*np.exp(-alpha*(x-xbeta)**2+1j*pbeta*(x-xbeta)/hbar)

xsorted = np.sort(x)
xmin = -3
xmax = 3
xminindex = np.where(xsorted>xmin)[0][0]
xmaxindex = np.where(xsorted>xmax)[0][0]

p = x.argsort()

fig = plt.figure(1)
plt.cla()
plt.ion()
plt.plot(xsorted[xminindex:xmaxindex],np.real(Phi[p])[xminindex:xmaxindex], 'b-')
plt.draw()

raw_input("")

dt = 1.0/float(N*10)
timesteps = 100

for i in range(timesteps):
    plt.cla()
    Phi = np.dot(U,np.dot(sclinalg.expm2(-1j*HDVR*dt/hbar),np.dot(Ucc,Phi)))
    plt.plot(xsorted[xminindex:xmaxindex],np.real(Phi[p])[xminindex:xmaxindex], 'b-')
    plt.draw()
    time.sleep(0.1)

plt.show()

