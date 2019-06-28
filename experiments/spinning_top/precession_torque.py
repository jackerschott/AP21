import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt, exp
import scipy.constants as cs
import os

import datproc.plot as dp

from damping import delta, d_delta, omega_fit

## General
output = __name__ == '__main__'

## Data
f_F = np.array([
  [680.0, 580.0, 375.0, 265.0],
  [675.0, 570.0, 370.0, 260.0],
  [680.0, 515.0, 415.0, 220.0],
  [640.0, 520.0, 370.0, 280.0]
])
d_f_F = np.array([
  [10.0, 10.0, 10.0, 10.0],
  [10.0, 10.0, 10.0, 10.0],
  [10.0, 10.0, 10.0, 10.0],
  [10.0, 10.0, 10.0, 10.0]
])
T_P = np.array([
  [45.4, 41.3, 29.2, 19.2],
  [41.1, 33.7, 23.9, 15.0],
  [29.2, 22.7, 17.6, 10.3],
  [21.8, 17.8, 13.2, 9.0]
])
d_T_P = np.array([
  [1.0, 1.0, 1.0, 1.0],
  [1.0, 1.0, 1.0, 1.0],
  [1.0, 1.0, 1.0, 1.0],
  [1.0, 1.0, 1.0, 1.0],
])

## Data processing
omega_F1 = 2 * pi * f_F
d_omega_F1 = 2 * pi * d_f_F

omega_F2 = omega_F1 * exp(-delta * T_P)
d_omega_F2 = omega_F2 * sqrt((d_omega_F1 / omega_F1)**2 + (d_delta * T_P)**2 + (d_T_P * delta)**2)

omega_F = 0.5 * (omega_F1 + omega_F2)
d_omega_F = sqrt(d_omega_F1**2 + d_omega_F2**2)

## Evaluation
if output:
  plt.subplots()
  plt.xlabel(r'$\omega_F$')
  plt.ylabel(r'$T_P$')

s, d_s = np.zeros(len(f_F)), np.zeros(len(f_F))
b, d_b = np.zeros(len(f_F)), np.zeros(len(f_F))
for i in range(len(f_F)):
  s[i], d_s[i], b[i], d_b[i] = dp.linreg(omega_F[i], T_P[i], d_T_P[i], d_omega_F[i])
  x_fit = dp.x_fit_like(omega_F[i])
  y_fit, y_u_fit = dp.linreg_lines(x_fit, s[i], d_s[i], b[i], d_b[i])

  dataPts, *_ = plt.errorbar(omega_F[i], T_P[i], d_T_P[i], d_omega_F[i], fmt='o')
  plt.plot(x_fit, y_fit, color=dataPts.get_color())
  plt.plot(x_fit, y_u_fit, color=dataPts.get_color(), ls='dashed')

if output:
  fig_folder_path = 'figures/spinning_top'
  if not os.path.exists(fig_folder_path):
    os.makedirs(fig_folder_path)

  fig_paths = dp.get_fig_paths(fig_folder_path, plt.get_fignums(), format='pdf')
  for fig_path in fig_paths:
    plt.savefig(fig_path, bbox_inches='tight', pad_inches=0.6)

  fig_paths = dp.get_fig_paths(fig_folder_path, plt.get_fignums(), format='pgf')
  for fig_path in fig_paths:
    plt.savefig(fig_path, bbox_inches='tight', pad_inches=0.6)

plt.show()
