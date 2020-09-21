import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

def FTCS(f, c, imax):
   fnew = np.zeros(imax+2, float)
   for i in range(1, imax+1):
      fnew[i] = f[i] - c*0.5*(f[i+1] - f[i-1])
   return fnew

def lax(f, c, imax):
   fnew = np.zeros(imax+2, float)
   for i in range(1, imax+1):
      fnew[i] = 0.5*(f[i+1] + f[i-1]) - c*0.5*(f[i+1] - f[i-1])
   return fnew

def upwind(f, c, imax):
   fnew = np.zeros(imax+2, float)
   for i in range(1, imax+1):
      fnew[i] = f[i] - c*(f[i] - f[i-1])
   return fnew

# for exercise eq. 11.37
def lax_wendroff(f, c, imax):
   fnew = np.zeros(imax+2, float)
   for i in range(1, imax+1):
      fnew[i] = f[i] - 0.5*c*(f[i+1] - f[i-1]) + c*c*0.5*(f[i+1] - 2*f[i] + f[i-1])
   return fnew

imax = 50
itrmax = 100

c = 0.5    # convection number
m = 5      # mode number (nuber of waves)
theta = 2*np.pi*m/imax   # phase difference in single grid interval

x = np.zeros(imax+2, float)     # initial condition
f = np.zeros(imax+2, float)
for i in range(0, imax+2):
   x[i] = i-1
   f[i] = np.cos((i-1)*theta)

xr = np.linspace(0, imax, 1000)   # exact solution for reference
fr = np.cos(theta*xr)

fig = plt.figure()
ims = []
plt.title('theta = %1.5f, c = %1.1f'%(theta,c))
plt.xlim(0,imax-1)
plt.xlabel('grid number')
plt.ylabel('f')
img = plt.plot(xr, fr, color='orange') \
+ plt.plot(x, f, marker='.', markersize=10, color='blue')
ims.append(img)

for itr in range(1, itrmax+1):
   f[0] = f[imax]      # periodic boundary condition
   f[imax+1] = f[1]    # periodic boundary condition
   f = lax_wendroff(f, c, imax)
   fr = np.cos(theta*(xr-c*itr))   # exact solution
   img = plt.plot(xr, fr, color='orange')\
   + plt.plot(x, f, marker='.', markersize=10, color='blue')
   ims.append(img)

ani = animation.ArtistAnimation(fig, ims, interval=50)
plt.show()
