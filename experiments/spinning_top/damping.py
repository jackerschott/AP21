import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt, exp, log
import os
import scipy.constants as cs
from scipy.optimize import curve_fit

import datproc.plot as dp
import datproc.print as dpr

## General
output = __name__ == '__main__'

## Data
f = np.array([650, 645, 583, 510, 472, 455, 432]) / cs.minute
d_f = np.array([5, 5, 5, 5, 5, 5, 5]) / cs.minute
t = np.array([0, 2, 4, 6, 8, 10, 12]) * cs.minute

## Data processing
omega = 2 * pi * f
d_omega = 2 * pi * d_f

## Evaluation
def omega_fit(t, delta, omega_0):
  return omega_0 * exp(-delta * t)

popt, pcov = curve_fit(omega_fit, t, omega, p0=(1/1000, 80), sigma=d_omega)
if output:
  plt.subplots(num=1)
  plt.xlabel(r'$t$ / s')
  plt.ylabel(r'$\omega$ / (1/s)')
  plt.yscale('log')
  lines, *_ = plt.errorbar(t, omega, d_omega, fmt='o')

  x_fit = dp.x_fit_like(t)
  y_fit = omega_fit(x_fit, *popt)
  plt.plot(x_fit, y_fit, color=lines.get_color())

delta, *_ = popt
d_delta = sqrt(pcov[0, 0])

t_h = log(2) / delta
d_t_h = t_h * d_delta / delta

if output:
  print(dpr.val(delta * cs.day, d_delta * cs.day, name='δ', unit='1/d', prefix=False, exp_to_fix=0))
  print(dpr.val(t_h / cs.minute, d_t_h / cs.minute, name='T_H', unit='min', prefix=False, exp_to_fix=0))

if output:
  fig_folder_path = 'figures/spinning_top'
  if not os.path.exists(fig_folder_path):
    os.makedirs(fig_folder_path)
  
  fig_paths = dp.get_fig_paths(fig_folder_path, plt.get_fignums(), format='pgf')
  for i, path in zip(plt.get_fignums(), fig_paths):
    plt.figure(i).savefig(path, bbox_inches='tight', pad_inches=0.0)
  
  fig_paths = dp.get_fig_paths(fig_folder_path, plt.get_fignums(), format='pdf')
  for i, path in zip(plt.get_fignums(), fig_paths):
    plt.figure(i).savefig(path, bbox_inches='tight', pad_inches=0.2)
