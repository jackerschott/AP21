import numpy as np
from numpy import sqrt
import scipy.constants as cs

import datproc.print as dpr

from cooling_unit import c_W, rho_W

## General
output = __name__ == '__main__'

## Data
I = 2.590 * 5
d_I = 0.010 * 5

U = 11.78
d_U = 0.03

J = 212.0 * cs.milli * cs.liter / cs.minute
d_J = 2.0 * cs.milli * cs.liter / cs.minute

f = np.array([353.0, 353.0, 353.0]) / cs.minute
d_f = np.array([1.0, 1.0, 1.0]) / cs.minute

T_out = 25.00
d_T_out = 0.10
T_in = 17.80
d_T_in = 0.10

A = np.array([15236, 15346, 15383]) * cs.hecto * cs.centi**3

## Data procession
f = np.mean(f)
d_f = sqrt(np.sum(d_f**2)) / len(d_f)

delta_T = T_out - T_in
d_delta_T = sqrt(d_T_in**2 + d_T_out**2)

d_A = np.std(A) / sqrt(len(A))
A = np.mean(A)

## Evaluation

Q_el = I * U / f
d_Q_el = Q_el * sqrt((d_I / I)**2 + (d_U / U)**2 + (d_f / f)**2)

Q_out = c_W * rho_W * delta_T * J / f
d_Q_out = Q_out * sqrt((d_delta_T / delta_T)**2 + (d_J / J)**2 + (d_f / f)**2)

W_pV = A
d_W_pV = d_A

Q_V = Q_el - Q_out - W_pV
d_Q_V = sqrt(d_Q_el**2 + d_Q_out**2 + d_W_pV**2)

P_el = Q_el * f
d_P_el = P_el * sqrt((d_Q_el / Q_el)**2 + (d_f / f)**2)

P_out = Q_out * f
d_P_out = P_out * sqrt((d_Q_out / Q_out)**2 + (d_f / f)**2)

P_pV = W_pV * f
d_P_pV = P_pV * sqrt((d_W_pV / W_pV)**2 + (d_f / f)**2)

eta_th = W_pV / Q_el
d_eta_th = eta_th * sqrt((d_W_pV / W_pV)**2 + (d_Q_el / Q_el)**2)

if output:
  print(dpr.val(Q_el, d_Q_el, name='Q_el', unit='J'))
  print(dpr.val(Q_out, d_Q_out, name='Q_out', unit='J'))
  print(dpr.val(W_pV, d_W_pV, name='W_pV', unit='J'))
  print(dpr.val(Q_V, d_Q_V, name='Q_V', unit='J'))
  print()
  print(dpr.val(P_el, d_P_el, name='P_el', unit='W'))
  print(dpr.val(P_out, d_P_out, name='P_out', unit='W'))
  print(dpr.val(P_pV, d_P_pV, name='P_pV', unit='W'))
  print()
  print(dpr.val(eta_th, d_eta_th, name='eta_th', exp_to_fix=0))
  print(dpr.val(eta_th * 100, d_eta_th * 100, name='eta_th', unit='%', exp_to_fix=0))
