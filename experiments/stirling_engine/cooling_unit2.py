import numpy as np
from numpy import sqrt
import scipy.constants as cs

import datproc.print as dpr

from cooling_unit import c_W, rho_W, UM, d_UM, IM, d_IM, eta, d_eta

from cooling_unit import f

## General
output = __name__ == '__main__'

## Data
lda_W = 335.0 / cs.gram

V = 1.0 * cs.milli * cs.liter

t1 = 335.0
d_t1 = 5.0
t2 = 565.0
d_t2 = 5.0

## Evaluation
delta_t = t2 - t1
d_delta_t = sqrt(d_t1**2 + d_t2**2)

Q2 = rho_W * V * lda_W

WM = UM * IM * delta_t
d_WM = Q2 * sqrt((d_UM / UM)**2 + (d_IM / IM)**2 + (d_delta_t / delta_t)**2)

eta_ = Q2 / WM
d_eta_ = eta_ * d_WM / WM

if output:
  print(dpr.val(delta_t, d_delta_t, name='Δt', unit='s'))
  print()
  print(dpr.val(Q2, name='Q2', unit='J'))
  print(dpr.val(WM, d_WM, name='WM', unit='J'))
  print()
  print(dpr.val(eta, d_eta, name='η', exp_to_fix=0))
  print(dpr.val(eta_, d_eta_, name='η\'', exp_to_fix=0))
  print(dpr.dev(eta, d_eta, eta_, d_eta_, name='η, η\''))
