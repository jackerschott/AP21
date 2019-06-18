import numpy as np

import datproc.print as dpr

## General
output = __name__ == '__main__'

T = np.array([37.0, 34.9, 36.1])
d_T = np.array([1.0, 1.0, 1.0])

if output:
  print(dpr.dev(T[0], d_T[0], T[1], d_T[1], name='T1, T2'))
  print(dpr.dev(T[0], d_T[0], T[2], d_T[2], name='T1, T3'))
  print(dpr.dev(T[1], d_T[1], T[2], d_T[2], name='T2, T3'))
