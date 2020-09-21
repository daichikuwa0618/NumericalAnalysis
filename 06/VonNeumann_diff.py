import numpy as np
import matplotlib.pyplot as plt

def FTCS(f, th, d):
   fm1 = f*np.exp(-1j*th)
   fp1 = f*np.exp( 1j*th)
   fnew = f + d*(fp1 - 2.0*f + fm1)
   return fnew

def BTCS(f, th, d):
   fm1 = f*np.exp(-1j*th)
   fp1 = f*np.exp( 1j*th)
   fold = -d*fp1 + (1.0+2.0*d)*f -d*fm1
   fnew = (f/fold)*f
   return fnew

d = 0.5    # diffusion number
imax = 100
theta = np.linspace(0, 2*np.pi, imax)
x = np.zeros(imax, float)
y = np.zeros(imax, float)

xc = np.cos(theta)   # unit circle for reference
yc = np.sin(theta)

for i in range(0, imax):
   f = 1.0 + 0j    # arbitrary complex number input
   fnew = FTCS(f, theta[i], d)
   G = fnew/f
   x[i] = G.real
   y[i] = G.imag

# plot setting
plt.xlabel('Re(G)')
plt.ylabel('Im(G)')
plt.plot(xc, yc, color='orange')
plt.plot(x, y, marker='.', markersize=10, color='blue')
plt.show()
