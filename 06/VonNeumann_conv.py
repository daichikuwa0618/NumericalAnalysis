import numpy as np
import matplotlib.pyplot as plt

def FTCS(f, th, c):
   fm1 = f*np.exp(-1j*th)
   fp1 = f*np.exp( 1j*th)
   fnew = f - c*0.5*(fp1 - fm1)
   return fnew

def lax(f, th, c):
   fm1 = f*np.exp(-1j*th)
   fp1 = f*np.exp( 1j*th)
   fnew = 0.5*(fp1 + fm1) - c*0.5*(fp1 - fm1)
   return fnew

def upwind(f, th, c):
   fm1 = f*np.exp(-1j*th)
   fnew = f - c*(f - fm1)
   return fnew

def lax_wendroff(f, th, c):
   fm1 = f*np.exp(-1j*th)
   fp1 = f*np.exp( 1j*th)
   fnew = f - c*0.5*(fp1-fm1) + c*c*0.5*(fp1-2*f+fm1)
   return fnew

c = 1.5    # convection number
imax = 100
theta = np.linspace(0, 2*np.pi, imax)
x = np.zeros(imax, float)
y = np.zeros(imax, float)

xc = np.cos(theta) # unit circle for reference
yc = np.sin(theta)

for i in range(0, imax):
   f = 1.0 + 0j    # arbitrary complex number input
   fnew = lax_wendroff(f, theta[i], c)
   G = fnew/f
   x[i] = G.real
   y[i] = G.imag

# plot setting
plt.xlabel('Re(G)')
plt.ylabel('Im(G)')
plt.plot(xc, yc, color='orange')
plt.plot(x, y, marker='.', markersize=10, color='blue')
plt.show()
