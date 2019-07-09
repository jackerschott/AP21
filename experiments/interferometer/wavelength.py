import numpy as np
from numpy import sqrt
import scipy.constants as cs

import datproc.print as dpr

## General
output = __name__ == '__main__'

## Data
lda_mfr = 532.0 * cs.nano
d_lda_mfr = 1.0 * cs.nano

s1 = np.array([0.0000, 0.3010, 0.6000, 0.9010, 1.2010]) * cs.milli
d_s1 = np.array([0.0010, 0.0010, 0.0010, 0.0010, 0.0010]) * cs.milli
s2 = np.array([3.0640, 3.3560, 3.6570, 3.9730, 4.2640]) * cs.milli
d_s2 = np.array([0.0010, 0.0010, 0.0010, 0.0010, 0.0010]) * cs.milli
d_s_sys = 9 * cs.micro

m = np.array([11170, 11142, 11163, 11187, 11149])
d_m = np.array([50, 50, 50, 50, 50])

## Evaluation
if output:
  print(dpr.val(np.std(m) / sqrt(len(s1)), name='Δm_stat'))
  print(dpr.val(sqrt(np.sum(d_s1**2 + d_s2**2)) / len(d_s1), name='Δs_read', unit='m'))

s = s2 - s1
d_s = np.std(s) / sqrt(len(s))
s = np.mean(s)

m = np.mean(m)
d_m = sqrt(np.sum(d_m**2)) / len(d_m)

if output:
  print(dpr.val(s, d_s, name='s', unit='m'))
  print(dpr.val(m, d_m, name='m'))
  print()

lda = 2 * s / m
d_lda = lda * sqrt(d_s**2 / s**2 + (d_m / m)**2)

if output:
  print(dpr.val(lda, d_lda, name='λ', unit='m'))
  print(dpr.dev(lda, d_lda, lda_mfr, d_lda_mfr, name='λ, λ_mfr'))
