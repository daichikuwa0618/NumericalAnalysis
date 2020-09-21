# Numerical Analysis calss Assignment
# created by Daichi Hayashi (B6TB1505) Oct. 07, 2019.
# Python script of "False Position method" (Regula Falsi)
import numpy as np
import matplotlib.pyplot as plt

def func(x):
  alpha = 40.0 # input angle [deg.]
  return 5.0/3.0*np.cos(np.deg2rad(alpha)) - 5.0/2.0*np.cos(np.deg2rad(x)) + 11.0/6.0 - np.cos(np.deg2rad(alpha-x))

def main():
  print('Start finding Root of function with the Regula Falsi...')
  n_max = 1000 # iteration max number
  num_a = 30.0 # smaller initial value [deg.]
  num_b = 40.0 # larger initial value [deg.]
  eps1  = 1e-9 # width threashold
  eps2  = 1e-9 # threashold

  if func(num_a)*func(num_b) > 0: # same sign -> initial value is wrong
    print("The initial values should be both side of root")
    return
  else:
    for n in range(n_max):
      fa    = func(num_a)
      fb    = func(num_b)
      num_c = (num_b*fa - num_a*fb)/(fa - fb)
      fc    = func(num_c)
      # print status
      print('loop {:>2d}, a = {:.7f}, f(a) = {:.9f}, b = {:.1f}, x = {:.7f}, f(x) = {:.9f}'.format(n+1,num_a,fa,num_b,fb,num_c,fc))
      if abs(fc) < eps2: # judge from function(c)
        print('Root has been found within accuracy (f(x) = {:5.1e}, loop = {:>2d}).'.format(eps2,n+1))
        return
      elif abs(num_b - num_a) < eps1: # judge from |b - a|, in almost case, this condition has NO meaning.
        print('Root has been found within the very small phi width (eps = {:5.1e}, iteration = {:>2d}).'.format(eps1,n+1))
        return
      if fa * fc > 0:
        num_a = num_c
      elif fb * fc > 0:
        num_b = num_c

if __name__ == '__main__':
  main()
