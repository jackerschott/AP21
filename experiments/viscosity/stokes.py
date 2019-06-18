import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt
import os
import scipy.constants as cs

import datproc.print as dpr
import datproc.plot as dp

output = __name__ == '__main__'
if output:
  print()

## Data
g = 9.80984
d_g = 0.00002

rho_peg = 1.1451 * cs.gram / cs.centi**3
d_rho_peg = 0.0004 * cs.gram / cs.centi**3
R = 75 * cs.milli / 2

r = np.array([ 4.5, 3.572, 4.0, 3.0, 2.5, 2.0, 1.5, 1.0, 0.75 ]) * cs.milli
d_r = 0.005 * r
rho_K = np.array([[1.360, 1.365], [1.375, 1.380], [1.355, 1.360], [1.375, 1.380], [1.375, 1.380],
                  [1.375, 1.380], [1.375, 1.380], [1.375, 1.380], [1.390, 1.395]]) * cs.gram / cs.centi**3

t = np.array([
  [13.3, 18.0, 16.5, 21.2, 29.0, 44.6, 56.3, 109.5, 110.4],
  [13.0, 17.8, 16.4, 21.5, 29.2, 43.8, 54.5, 108.8, 104.7],
  [12.8, 18.0, 16.4, 21.0, 29.1, 44.1, 56.5, 111.2, 113.6],
  [12.9, 17.8, 16.0, 21.4, 29.1, 43.5, 56.9, 109.7, 118.6],
  [13.0, 17.9, 16.3, 21.1, 29.1, 44.2, 56.9, 111.2, 113.6]
])
d_t = np.array([
  [0.3, 0.3, 0.3, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0],
  [0.3, 0.3, 0.3, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0],
  [0.3, 0.3, 0.3, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0],
  [0.3, 0.3, 0.3, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0],
  [0.3, 0.3, 0.3, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0]
])
s = np.array([450, 450, 450, 400, 400, 400, 300, 300, 200]) * cs.milli

T = 24.45
d_T = 0.05

## Data processing
t = np.mean(t, axis=0)
d_t = sqrt(np.sum(d_t**2, axis=0)) / np.size(d_t, axis=0)

d_rho_K = 0.5 * np.abs(np.diff(rho_K, axis=1)[:,0])
rho_K = 0.5 * np.sum(rho_K, axis=1)

## Evaluation
r2 = r**2
d_r2 = 2 * r * d_r
v = s / t
d_v = v * d_t / t

lda = 1 + 2.1 * r / R
d_lda = 2.0 * d_r / R
v_c = v * lda
d_v_c = sqrt(d_v**2 + d_lda**2)
if output:
  print(dpr.tbl([
    dpr.lst(r, d_r, name='r', unit='m'),
    dpr.lst(lda, d_lda, name='λ'),
    dpr.lst(v, d_v, name='v', unit='m / s'),
    dpr.lst(v_c, d_v_c, name='v_corr', unit='m / s'),
    dpr.dev(v, d_v, v_c, d_v_c, name='v, v_c')
  ]))
  print()

v_rho_ratio = v / (rho_K - rho_peg)
d_v_rho_ratio = v_rho_ratio * sqrt((d_v / v)**2 + (d_rho_K**2 + d_rho_peg**2) / (rho_K - rho_peg)**2)

if output:
  plt.subplots(num=1)
  plt.xlabel(r'$r^2$ / mm$^2$')
  plt.ylabel(r'$\frac{v}{\varrho_K - \varrho_\textnormal{PEG}} / \frac{\textnormal{cm}^4}{\textnormal{g} \, \textnormal{s}}$')
  plt.errorbar(*dp.to_units(r2, v_rho_ratio, d_v_rho_ratio, d_r2, x_unit=cs.milli**2, y_unit=(cs.centi**4 / cs.gram)), fmt='o')

v[:7] = v_c[:7]
d_v[:7] = d_v_c[:7]
v_rho_ratio = v / (rho_K - rho_peg)
d_v_rho_ratio = v_rho_ratio * sqrt((d_v / v)**2 + (d_rho_K**2 + d_rho_peg**2) / (rho_K - rho_peg)**2)
slope1, d_slope1, itc1, d_itc1 = dp.linreg(r2, v_rho_ratio, d_v_rho_ratio, d_r2)

if output:
  lines, *_ = plt.errorbar(*dp.to_units(r2, v_rho_ratio, d_v_rho_ratio, d_r2, x_unit=cs.milli**2, y_unit=(cs.centi**4 / cs.gram)), fmt='o')

  x_line = dp.x_fit_like(r2)
  y_line, y_uline = dp.linreg_lines(x_line, slope1, d_slope1, itc1, d_itc1)
  plt.plot(*dp.to_units(x_line, y_line, x_unit=(cs.milli**2), y_unit=(cs.centi**4 / cs.gram)),
    label='Fit', color=lines.get_color())
  plt.plot(*dp.to_units(x_line, y_uline, x_unit=(cs.milli**2), y_unit=(cs.centi**4 / cs.gram)),
    label='Fit uncertainty', color=lines.get_color(), ls='dashed')
  plt.legend()

if output:
  print(dpr.val(slope1 / (cs.centi**2 / cs.gram), d_slope1 / (cs.centi**2 / cs.gram),
    name='slope1', unit='cm^2 / g', prefix=False))
  print()

eta = 2 / 9 * g / slope1
d_eta = eta * sqrt((d_g / g)**2 + (d_slope1 / slope1)**2)

if output:
  print(dpr.val(eta, d_eta, name='η', unit='Pa s'))

v_theo = 2/9 * g * (rho_K - rho_peg) * r2 / eta
d_v_theo = v_theo * sqrt((d_g / g)**2 + (d_rho_K**2 + d_rho_peg**2) / (rho_K - rho_peg)**2 + (d_r2 / r2)**2 + (d_eta / eta)**2)
v__v_theo = v / v_theo
d_v__v_theo = v__v_theo * sqrt((d_v / v)**2 + (d_v_theo / v_theo)**2)

Re = rho_peg * v * R / eta
d_Re = Re * sqrt((d_rho_peg / rho_peg)**2 + (d_v / v)**2 + (d_eta / eta)**2)

if output:
  plt.subplots(num=2)
  plt.xlabel(r'$v / v_\textnormal{theo}$')
  plt.ylabel(r'Re')
  plt.yscale('log')
  plt.errorbar(v__v_theo, Re, d_Re, d_v__v_theo, fmt='o')

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
