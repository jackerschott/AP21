import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt
import os
import scipy.constants as cs

import datproc.plot as dpl
import datproc.print as dpr

from precession_I_z import I_z, d_I_z

## General
output = __name__ == '__main__'

## Data
f_F = np.array([
  775.0, 625.0, 420.0, 500.0, 800.0,
  900.0, 1030.0, 1115.0, 1265.0, 530.0
]) / cs.minute
d_f_F = np.array([
  10.0, 10.0, 10.0, 10.0, 10.0,
  10.0, 10.0, 10.0, 10.0, 10.0
]) / cs.minute
f_N = np.array([
  750.0, 610.0, 395.0, 490.0, 790.0,
  885.0, 955.0, 1055.0, 1175.0, 515.0
]) / cs.minute
d_f_N = np.array([
  10.0, 10.0, 10.0, 10.0, 10.0,
  10.0, 10.0, 10.0, 10.0, 10.0
]) / cs.minute

## Data procesing
omega_F = 2 * pi * f_F
d_omega_F = 2 * pi * d_f_F
omega_N = 2 * pi * f_N
d_omega_N = 2 * pi * d_f_N

## Evaluation
if output:
  plt.subplots(num=4)
  plt.xlabel(r'$\omega_N$ / (1/s)')
  plt.ylabel(r'$\omega_F$ / (1/s)')

s, d_s, b, d_b = dpl.linreg(omega_N, omega_F, d_omega_F, d_omega_N)
if output:
  print(dpr.val(s, d_s, name='s'))
if output:
  x_fit = dpl.x_fit_like(omega_N)
  y_fit, y_u_fit = dpl.linreg_lines(x_fit, s, d_s, b, d_b)

  dataPts, *_ = plt.errorbar(omega_N, omega_F, d_omega_F, d_omega_F, fmt='o')
  plt.plot(x_fit, y_fit, label='Fit', color=dataPts.get_color())
  plt.plot(x_fit, y_u_fit, label='Fit uncertainty', color=dataPts.get_color(), ls='dashed')
  plt.legend()

I_x = I_z * s
d_I_x = I_x * sqrt((d_I_z / I_z)**2 + (d_s / s)**2)

if output:
  print(dpr.val(I_x / (cs.gram * cs.centi**2), d_I_x / (cs.gram * cs.centi**2),
    name='I_x', unit='g cm^2'))

if output:
  fig_folder_path = 'figures/spinning_top'
  if not os.path.exists(fig_folder_path):
    os.makedirs(fig_folder_path)

  fig_paths = dpl.get_fig_paths(fig_folder_path, plt.get_fignums(), format='pdf')
  for fig_path in fig_paths:
    plt.savefig(fig_path, bbox_inches='tight', pad_inches=0.6)

  fig_paths = dpl.get_fig_paths(fig_folder_path, plt.get_fignums(), format='pgf')
  for fig_path in fig_paths:
    plt.savefig(fig_path, bbox_inches='tight', pad_inches=0.0)

plt.show()
