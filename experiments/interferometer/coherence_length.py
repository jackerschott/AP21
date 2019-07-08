import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt, exp
import os
import scipy.constants as cs
from scipy.optimize import curve_fit
from scipy.signal import find_peaks_cwt

import datproc.print as dpr
import datproc.plot as dpl

## General
output = __name__ == '__main__'

def gauss(t, A, mu, sigma):
  return A * exp(-0.5 * (t - mu)**2 / sigma**2) / sqrt(2 * pi * sigma**2)

## Data
lda = 532.0 * cs.nano
d_lda = 1.0 * cs.nano
v_mov = 0.1 * cs.milli

data = np.genfromtxt('data/interferometer/signal.csv', delimiter=',', skip_header=18)
t = data[:,3]
U = data[:,4]

## Evaluation
t_range = [-0.035, 0.065]
if output:
  plt.subplots(num=4)
  plt.xlabel(r'$t$ / s')
  plt.ylabel(r'$U$ / V')
  plt.xlim(*t_range)
  plt.plot(t, U)

i_U_max = find_peaks_cwt(U, range(1,30), noise_perc=20)
i_U_max = i_U_max[U[i_U_max] > 0.0]
i_U_max = i_U_max[t[i_U_max] > t_range[0]]
i_U_max = i_U_max[t[i_U_max] < t_range[1]]

popt, pcov = curve_fit(gauss, t[i_U_max], U[i_U_max], p0=(2.5e-3, 0.015, 0.015))
d_popt = sqrt(np.diag(pcov))
if output:
  x_fit = dpl.x_fit_like(t[i_U_max])
  y_fit = gauss(x_fit, *popt)

  data_pts, *_ = plt.plot(t[i_U_max], U[i_U_max], marker='o', ls='None')
  plt.plot(x_fit, y_fit, color=data_pts.get_color())

if output:
  print(dpr.val(popt[2], d_popt[2], name='Δt', unit='s'))

delta_s = v_mov * popt[2]
d_delta_s = v_mov * d_popt[2]

delta_lda = lda**2 / (2 * pi * delta_s)
d_delta_lda = delta_lda * sqrt((2 * d_lda / lda)**2 + (d_delta_s / delta_s)**2)

if output:
  print(dpr.val(delta_s, d_delta_s, name='Δs', unit='m'))
  print(dpr.val(delta_lda, d_delta_lda, name='Δλ', unit='m'))

L = lda**2 / delta_lda
d_L = L * sqrt((2 * d_lda / lda)**2 + (d_delta_lda / delta_lda)**2)

if output:
  print(dpr.val(L, d_L, name='L', unit='m'))

if output:
  fig_folder_path = 'figures/interferometer'
  if not os.path.exists(fig_folder_path):
    os.makedirs(fig_folder_path)
  for i in plt.get_fignums():
    plt.figure(i).savefig(os.path.join(fig_folder_path, 'fig' + str(i) + '.pgf'), bbox_inches='tight', pad_inches=0.0)
    plt.figure(i).savefig(os.path.join(fig_folder_path, 'fig' + str(i) + '.pdf'))
  plt.show()
