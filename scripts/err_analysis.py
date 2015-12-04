import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

A = np.loadtxt('prec.steps_keys', skiprows=1)
B = np.vstack([range(A.size), A]).T

#Check if error is indeed in the same form as previously
def func(x, a, b, c):
  return a*(b/(b+x))**c

x, y = zip(*B)

beg = 5
wb = 1e9
end = 10
we = 1e-3

sigma = [wb]*beg+[1]*(len(y) - (beg+end)) + [we]*end

popt, pcov = curve_fit(func, x, y, p0=[1., 7., 2.], sigma=sigma)

print popt
print pcov

plt.loglog(x, y)
plt.loglog(x, [func(v, *popt) for v in x])
plt.title('Fit : 2*(4.5/(4.5+x))'.format(*popt))
plt.ylabel('Error (m^3)')
plt.xlabel('Iterations')
#plt.loglog(x, [2*(6/(6+v)) for v in x])
plt.show()
