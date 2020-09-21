import numpy as np

def sor(a,b,eps,imax):
  n = b.size # dimension of equation
  x = np.zeros(n, float) # zero matrix

  for it in range(0,imax):
    dxmax = 0.0
    omega = 1.1 # acceleration coefficient
    for i in range(0,n):
      res = b[i] # residual
      for j in range(0,n):
        res -= a[i,j]*x[j]
      dxmax = max(abs(res),dxmax)
      x[i] = x[i] + omega*res/a[i,i]
    print(it, dxmax)
    if dxmax < eps: return x

def main():
  a = np.array([[ 4.0,-1.0, 0.0, 1.0, 0.0],
                [-1.0, 4.0,-1.0, 0.0, 1.0],
                [ 0.0,-1.0, 4.0,-1.0, 0.0],
                [ 1.0, 0.0,-1.0, 4.0,-1.0],
                [ 0.0, 1.0, 0.0,-1.0, 4.0]])
  b = np.array([100.0,100.0,100.0,100.0,100.0])
  imax = 100
  eps = 1e-13

  print("A, b = ")
  print(a,",",b)
  sol = sor(a,b,eps,imax)
  print(sol)

if __name__=='__main__':
  main()
