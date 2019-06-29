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
f = np.array([
  600.0, 480.0, 410.0, 350.0, 445.0,
  375.0, 300.0, 530.0, 505.0, 570.0,
])
d_f = np.array([
  10.0, 10.0, 10.0, 10.0, 10.0,
  10.0, 10.0, 10.0, 10.0, 10.0
])
t = np.array([
  31.6, 35.9, 41.1, 50.7, 40.3,
  45.5, 54.6, 32.8, 35.1, 32.4
])
d_t = np.array([
  2.0, 2.0, 2.0, 2.0, 2.0,
  2.0, 2.0, 2.0, 2.0, 2.0
])

## Data processing
omega_F = 2 * pi * f
d_omega_F = omega_F * d_f / f
Omega = 2 * pi / t
d_Omega = Omega * d_t / t

## Evaluation
s, d_s, b, d_b = dpl.linreg(Omega, omega_F, d_omega_F, d_Omega)
if output:
  x_fit = dpl.x_fit_like(Omega)
  y_fit, y_u_fit = dpl.linreg_lines(x_fit, s, d_s, b, d_b)

  dataPts, *_ = plt.errorbar(Omega, omega_F, d_omega_F, d_Omega, fmt='o')
  plt.plot(x_fit, y_fit, color=dataPts.get_color(), label='Fit')
  plt.plot(x_fit, y_u_fit, color=dataPts.get_color(), label='Fit uncertainty', ls='dashed')
  plt.legend()

delta_I = I_z / (s - 1)
d_delta_I = delta_I * sqrt((d_I_z / I_z)**2 + (d_s / (s - 1))**2)

I_x = delta_I + I_z
d_I_x = sqrt(d_delta_I**2 + d_I_z**2)

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