import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt
import os
import scipy.constants as cs

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
slope1, d_slope1, itc1, d_itc1 = dp.linreg(V, t, d_t)
if output:
  plt.subplots(num=3)
  plt.xlabel(r'V / cm$^3$')
  plt.ylabel(r't / s')
  lines, *_ = plt.errorbar(*dp.to_units(V, t, d_t, x_unit=cs.centi**3), fmt='o')

  x_line = dp.x_fit_like(V)
  y_line, y_uline = dp.linreg_lines(x_line, slope1, d_slope1, itc1, d_itc1)
  plt.plot(*dp.to_units(x_line, y_line, x_unit=cs.centi**3), label='Fit', color=lines.get_color())
  plt.plot(*dp.to_units(x_line, y_uline, x_unit=cs.centi**3), label='Fit uncertainty', color=lines.get_color(), ls='dashed')

if output:
  print(dpr.val(slope1 * cs.centi**3, d_slope1 * cs.centi**3, name='slope1', unit='s / cm^3'))

J = 1 / slope1
d_J = J * d_slope1 / slope1

if output:
  print(dpr.val(J / cs.milli**3, d_J / cs.milli**3, name='J', unit='mm^3 / s'))

p_tube = h_mean * rho_peg * g
d_p_tube = p_tube * sqrt((d_h_mean / h_mean)**2 + (d_rho_peg / rho_peg)**2 + (d_g / g)**2)

if output:
  print(dpr.val(p_tube, d_p_tube, name='p_tube', unit='Pa'))

eta = pi * p_tube * R**4 / (8 * J * L)
d_eta = eta * sqrt((d_p_tube / p_tube)**2 + (4 * d_R / R)**2 + (d_J / J)**2 + (d_L / L)**2)

if output:
  print(dpr.val(eta, d_eta, name='η', unit='Pa s'))

if output:
  fig_folder_path = 'figures/viscosity'
  if not os.path.exists(fig_folder_path):
    os.makedirs(fig_folder_path)

  fig_paths = dp.get_fig_paths(fig_folder_path, plt.get_fignums(), format='pgf')
  for i, path in zip(plt.get_fignums(), fig_paths):
    plt.figure(i).savefig(path, bbox_inches='tight', pad_inches=0.0)
  
  fig_paths = dp.get_fig_paths(fig_folder_path, plt.get_fignums(), format='pdf')
  for i, path in zip(plt.get_fignums(), fig_paths):
    plt.figure(i).savefig(path, bbox_inches='tight', pad_inches=0.2)
