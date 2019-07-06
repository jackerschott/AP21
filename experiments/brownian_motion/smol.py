import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt, exp
import os
import scipy.constants as cs
from scipy.optimize import curve_fit
from scipy.stats import norm

import datproc.plot as dpl
import datproc.print as dpr

import stokes as st
from stokes import a, d_a, T, d_T, eta, d_eta, t, r2_

## General
output = __name__ == '__main__'

## Evaluation
r2_cum = np.cumsum(r2_)
t_ = np.append(t[:67], t[68:-1] - (t[67] - t[66]))
popt, pcov = curve_fit(lambda x, s, b: s * x + b, t_, r2_cum)
d_popt = sqrt(np.diag(pcov))
if output:
  x_fit = dpl.x_fit_like(t_)
  y_fit, y_u_fit = dpl.linreg_lines(x_fit, popt[0], d_popt[0], popt[1], d_popt[1])

  plt.subplots(num=3)
  plt.xlabel(r'$t$ / s')
  plt.ylabel(r'$r^2$ / $\mu$m$^2$')
  plt.plot(t_, r2_cum / cs.micro**2, marker='o', ls='None')
  plt.plot(x_fit, y_fit / cs.micro**2)

if output:
  print(dpr.val(popt[0] / cs.micro**2, d_popt[0] / cs.micro**2, name='slope', unit='μm² / s'))
  print()

D = popt[0] / 4
d_D = d_popt[0] / 4

k = 6 * pi * eta * a * D / T
d_k = k * sqrt((d_eta / eta)**2 + (d_a / a)**2 + (d_D / D)**2 + (d_T / T)**2)

if output:
  print(dpr.val(D / cs.micro**2, d_D / cs.micro**2, name='D', unit='μm²', prefix=False, exp_to_fix=0))
  print(dpr.val(k, d_k, name='k', unit='J / K', prefix=False))
  print(dpr.dev(k, d_k, cs.k, name='k, k_lit'))
  print()
  print(dpr.dev(D, d_D, st.D, st.d_D, name='D, D_st'))
  print(dpr.dev(k, d_k, st.k, st.d_k, name='k, k_st'))

## Save plots
if output:
  fig_folder_path = 'figures/brownian_motion'
  if not os.path.exists(fig_folder_path):
    os.makedirs(fig_folder_path)

  for i in plt.get_fignums():
    plt.figure(i).savefig(os.path.join(fig_folder_path, 'fig%i.pgf' % i), bbox_inches='tight', pad_inches=0.0)
    plt.figure(i).savefig(os.path.join(fig_folder_path, 'fig%i.pdf' % i), bbox_inches='tight')
