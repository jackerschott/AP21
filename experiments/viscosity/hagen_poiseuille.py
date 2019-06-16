import numpy as np
import scipy.constants as cs
from numpy import pi, sqrt

import datproc.print as dpr
import datproc.plot as dp

from stokes import rho_peg, d_rho_peg, g, d_g

output = __name__ == '__main__'
if output:
  print()

## Data
R = 1.5 * cs.milli / 2
d_R = 0.01 * cs.milli / 2
L = 100 * cs.milli
d_L = 0.5 * cs.milli

hI = 548.0 * cs.milli
d_hI = 1.0 * cs.milli
hF = 542.0 * cs.milli
d_hF = 1.0 * cs.milli

V = np.array([5.0, 10.0, 15.0, 20.0, 25.0]) * cs.centi**3
t = np.array([120.0, 240.0, 360.0, 480.0, 600.0])
d_t = np.array([5.0, 5.0, 5.0, 5.0, 5.0])

T = 25.0
d_T = 0.05

## Data processing
h_mean = 0.5 * (hI + hF)
d_h_mean = 0.5 * (hI - hF)

## Evaluation
dp.initplot(num=3, xlabel='V / cm$^3$', ylabel='t / s')
slope1, d_slope1, itc1, d_itc1 = dp.linreg(V, t, d_t, x_unit=(cs.centi**3), plot=True)

if output:
  print(dpr.val('slope1', slope1 * cs.centi**3, d_slope1 * cs.centi**3, unit='s / cm^3'))

J = 1 / slope1
d_J = J * d_slope1 / slope1

if output:
  print(dpr.val('J', J / cs.milli**3, d_J / cs.milli**3, unit='mm^3 / s'))

p_tube = h_mean * rho_peg * g
d_p_tube = p_tube * sqrt((d_h_mean / h_mean)**2 + (d_rho_peg / rho_peg)**2 + (d_g / g)**2)

if output:
  print(dpr.val('p_tube', p_tube, d_p_tube, unit='Pa'))

eta = pi * p_tube * R**4 / (8 * J * L)
d_eta = eta * sqrt((d_p_tube / p_tube)**2 + (4 * d_R / R)**2 + (d_J / J)**2 + (d_L / L)**2)

if output:
  print(dpr.val('η', eta, d_eta, unit='Pa s'))

  dp.savefigs('figures/viscosity', format='pgf')
  dp.savefigs('figures/viscosity', format='pdf')
