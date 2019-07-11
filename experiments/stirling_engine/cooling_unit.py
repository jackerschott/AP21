import numpy as np
from numpy import sqrt
import scipy.constants as cs

import datproc.print as dpr

## General
output = __name__ == '__main__'

## Data
c_W = 4180
rho_W = 998.2

IH = 1.06
d_IH = 0.01

UH = 5.35
d_UH = 0.01

IM = 1.70
d_IM = 0.05

UM = 23.90
d_UM = 0.10

J = 214.0 * cs.milli * cs.liter / cs.minute
d_J = 2.0 * cs.milli * cs.liter / cs.minute

delta_T = 3.25
d_delta_T = 0.07

f = 285.0 / cs.minute
d_f = 2.0 / cs.minute

## Evaluation
Q1 = c_W * rho_W * delta_T * J / f
d_Q1 = Q1 * sqrt((d_delta_T / delta_T)**2 + (d_J / J)**2 + (d_f / f)**2)

Q2 = UH * IH / f
d_Q2 = Q2 * sqrt((d_UH / UH)**2 + (d_IH / IH)**2 + (d_f / f)**2)

WM = UM * IM / f
d_WM = WM * sqrt((d_UM / UM)**2 + (d_IM / IM)**2 + (d_f / f)**2)

eta = Q2 / WM
d_eta = eta * sqrt((d_Q2 / Q2)**2 + (d_WM / WM)**2)

Q1_ = Q2 + WM
d_Q1_ = sqrt(d_Q2**2 + d_WM**2)

if output:
  print(dpr.val(Q1, d_Q1, name='Q1', unit='J'))
  print(dpr.val(Q2, d_Q2, name='Q2', unit='J'))
  print(dpr.val(WM, d_WM, name='WM', unit='J'))
  print()
  print(dpr.val(eta, d_eta, name='η', exp_to_fix=0))
  print(dpr.val(eta * 100, d_eta * 100, name='eta_th', unit='%', exp_to_fix=0))
  print()
  print(dpr.val(Q1, d_Q1, name='Q1', unit='J'))
  print(dpr.val(Q1_, d_Q1_, name='Q2 + WM', unit='J'))
  print(dpr.dev(Q1, d_Q1, Q1_, d_Q1_, name='Q1, Q2 + WM'))
