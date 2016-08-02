import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, leastsq, minimize
from copy import copy

A = np.loadtxt('../logs/prec.steps_keys', skiprows=1)
B = np.vstack([range(A.size), A]).T

#Check if error is indeed in the same form as previously
def func(x, a, b, c, d, e, f):
  return f+a/(b+(c+d*x)**e)

def slsqp_func(params, x, y):
  points = func(x, params[0], params[1], params[2], params[3], params[4], params[5])
  return np.linalg.norm(points - y)

def slsqp_func_pen(params, x, y):
  points = func(x, params[0], params[1], params[2], params[3], params[4], params[5])
  try:
    neg_points = points[points < y]
  except RuntimeWarning:
    neg_points = np.zeros(y.shape)
  return np.linalg.norm(points - y) + 0.001*neg_points.size

def residuals(p, x, y):
  points = func(x, *p)
  neg_points = copy(points)
  neg_points[points < y] = 1
  neg_points[points >= y] = 0
  #neg_points = neg_points*0
  #neg_points = np.zeros(y.shape)
  #print y.shape, points.shape, neg_points.shape
  return y - points - neg_points

x, y = B[:, 0], B[:, 1]

beg = 5
wb = 1e9
end = 10
we = 1e-3

sigma = [wb]*beg+[1]*(len(y) - (beg+end)) + [we]*end

popt, pcov = curve_fit(func, x, y, p0=[1., 7., 2., 2., 1., 0.])
#popt1, pcov1 = leastsq(func=residuals, args=(x, y), x0=[1., 7., 2., 2., 1., 0.])

cons = ({'type': 'ineq',
         'fun': lambda p: func(x, *p) - y+0.8})

res = minimize(slsqp_func, [1., 7., 2., 2., 1., 0.], args=(x, y),
               method='SLSQP', options={'disp': True}, constraints=cons)

print pcov
print
print popt

print res

plt.loglog(x, y, label='data')
plt.loglog(x, func(x, *popt), 'r--, ', label='fit unconstrained')
plt.loglog(x, func(x, *res['x']), 'g', label='fit constrained')
plt.title('Fit')
plt.ylabel('Error (m^3)')
plt.xlabel('Iterations')
plt.legend()
#plt.loglog(x, [2*(6/(6+v)) for v in x])
plt.show()
