import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt
import os
import scipy.constants as cs

import datproc.plot as dpl
import datproc.print as dpr

## General
output = __name__ == '__main__'

## Data
T0 = cs.zero_Celsius
p0 = cs.atm

n0_lit = 1.00028

lda = 532.0 * cs.nano
d_lda = 1.0 * cs.nano

a = 50.00 * cs.milli
d_a = 0.05 * cs.milli

T = 24.8 + cs.zero_Celsius
d_T = 0.3

n = np.array([0, 5, 10, 15, 20])
p_mano = np.array([
  [735.0, 680.0, 605.0, 535.0, 460.0],
  [730.0, 670.0, 595.0, 525.0, 450.0],
  [730.0, 670.0, 600.0, 525.0, 455.0]
]) * cs.torr

## Data procession
n = n[1:]
p = p_mano[:,0:1] - p_mano[:,1:]
d_p = p / (5 * cs.torr) * 0.006

popts = np.zeros((len(p),2))
d_popts = np.zeros((len(p),2))
for i in range(len(p)):
  popt = dpl.linreg(n, p[i], d_p[i])
  popts[i] = popt[0], popt[2]
  d_popts[i] = popt[1], popt[3]

if output:
  plt.subplots(num=1)
  for i in range(len(popts)):
    x_fit = dpl.x_fit_like(n)
    y_fit, y_u_fit = dpl.linreg_lines(x_fit,
      popts[i][0], d_popts[i][0], popts[i][1], d_popts[i][1])
    
    dataPts, *_ = plt.errorbar(n, p[i], d_p[i], fmt='o')
    plt.plot(x_fit, y_fit, color=dataPts.get_color())
    plt.plot(x_fit, y_u_fit, color=dataPts.get_color(), ls='dashed')

if output:
  print(dpr.tbl([
    dpr.lst(popts[:,0], d_popts[:,0], name='slope', unit='Pa')
  ]))

n0 = p0 * T * lda / (2 * T0 * a * popts[:,0]) + 1.0
d_n0 = n0 * sqrt((d_T / T)**2 + (d_lda / lda)**2 + (d_a / a)**2 + (d_popts[:,0] / popts[:,0])**2)

if output:
  print(dpr.tbl([
    dpr.lst(n0, d_n0, name='n0', prefix=False, exp_to_fix=0)
  ]))

n0 = np.mean(n0)
d_n0 = sqrt(np.sum(d_n0**2)) / len(d_n0)

if output:
  print(dpr.val(n0, d_n0, name='n0', prefix=False, exp_to_fix=0))
  print(dpr.dev(n0, d_n0, n0_lit, name='n0, n0_lit'))

if output:
  fig_folder_path = 'figures/interferometer'
  if not os.path.exists(fig_folder_path):
    os.makedirs(fig_folder_path)
  for i in plt.get_fignums():
    plt.figure(i).savefig(os.path.join(fig_folder_path, 'fig' + str(i) + '.pgf'), bbox_inches='tight', pad_inches=0.0)
    plt.figure(i).savefig(os.path.join(fig_folder_path, 'fig' + str(i) + '.pdf'))
  plt.show()
