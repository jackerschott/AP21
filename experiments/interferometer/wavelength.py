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

n = np.array([11170, 11142, 11163, 11187, 11149])
d_n = np.array([50, 50, 50, 50, 50])

## Evaluation
s = s2 - s1
d_s = sqrt((np.std(s) / sqrt(len(s)))**2 + d_s_sys**2)
s = np.mean(s)

lda = 2 * s / n
d_lda = lda * sqrt(d_s**2 / s**2 + (d_n / n)**2)

lda = np.mean(lda)
d_lda = sqrt(np.sum(d_lda**2)) / len(d_lda)

if output:
  print(dpr.val(lda, d_lda, name='λ', unit='m'))
  print(dpr.dev(lda, d_lda, lda_mfr, d_lda_mfr, name='λ, λ_mfr'))
