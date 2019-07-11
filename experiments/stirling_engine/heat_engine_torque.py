import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt
import os
import scipy.constants as cs

import datproc.print as dpr

from cooling_unit import c_W, rho_W

## General
output = __name__ == '__main__'

## Data
l = 25.0 * cs.centi

I = np.array([2.610, 2.600, 2.580, 2.560])
d_I = np.array([0.010, 0.010, 0.010, 0.010])

U = np.array([11.76, 11.76, 11.68, 11.64])
d_U = np.array([0.03, 0.03, 0.03, 0.03])

f = np.array([
  [222.0, 227.0, 223.0],
  [266.0, 277.0, 261.0],
  [278.0, 283.0, 275.0],
  [293.0, 294.5, 325.0]
]) / cs.minute
d_f = np.array([
  [2.0, 2.0, 2.0],
  [2.0, 2.0, 2.0],
  [2.0, 2.0, 2.0],
  [2.0, 2.0, 2.0]
]) / cs.minute

A = np.array([
  [24351, 24087, 24003],
  [22215, 21896, 21890],
  [20112, 19533, 19892],
  [17099, 16365, 16847]
]) * cs.hecto * cs.centi**3

F = np.array([0.93, 0.65, 0.48, 0.25])
d_F = np.array([0.15, 0.10, 0.10, 0.15])

## Data procession
f = np.mean(f, axis=1)
d_f = sqrt(np.sum(d_f**2, axis=1)) / np.size(d_f, axis=1)

d_A = np.std(A, axis=1) / sqrt(np.size(A, axis=1))
A = np.mean(A, axis=1)

## Evaluation
Q_el = I * U / f
d_Q_el = Q_el * sqrt((d_I / I)**2 + (d_U / U)**2 + (d_f / f)**2)

W_pV = A
d_W_pV = d_A

W_D = pi * l * F
d_W_D = W_D * d_F / F

eta_th = W_pV / Q_el
d_eta_th = eta_th * sqrt((d_W_pV / W_pV)**2 + (d_Q_el / Q_el)**2)

eta_eff = W_D / Q_el
d_eta_eff = eta_eff * sqrt((d_W_D / W_D)**2 + (d_Q_el / Q_el)**2)

if output:
  plt.subplots(num=1)
  plt.xlabel(r'$f$ / Hz')
  plt.ylabel(r'$\eta_\textnormal{th}$, $\eta_\textnormal{eff}$')
  plt.errorbar(f, eta_th, d_eta_th, d_f, fmt='o', ls='-')
  plt.errorbar(f, eta_eff, d_eta_eff, d_f, fmt='o', ls='-')

if output:
  fig_folder_path = 'figures/stirling_engine'
  if not os.path.exists(fig_folder_path):
    os.makedirs(fig_folder_path)
  for i in plt.get_fignums():
    plt.figure(i).savefig(os.path.join(fig_folder_path, 'fig%i.pgf' % i), bbox_inches='tight', pad_inches=0.0)
    plt.figure(i).savefig(os.path.join(fig_folder_path, 'fig%i.pdf' % i))
  plt.show()
