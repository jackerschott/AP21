import numpy as np
import scipy.constants as cs
from numpy import pi, sqrt

import datproc.print as dpr
import general as gen


## Data
tl = np.array([[1.23, 14.81], [1.26, 15.81], [1.45, 16.82]])
tr = np.array([[1.92, 15.49], [1.96, 16.53], [2.20, 17.58]])
d_tl = np.array([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
d_tr = np.array([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
f = np.array([0.736, 0.686, 0.649])
d_f = np.array([0.009, 0.01, 0.009])


## Data preparation
Tl = np.zeros(3)
d_Tl = np.zeros(3)
Tr = np.zeros(3)
d_Tr = np.zeros(3)
for i in range(3):
  Tl[i] = (tl[i][1] - tl[i][0]) / 10.0
  d_Tl[i] = sqrt(d_tl[i][1]**2 + d_tl[i][0]**2) / 10.0
  Tr[i] = (tr[i][1] - tr[i][0]) / 10.0
  d_Tr[i] = sqrt(d_tr[i][1]**2 + d_tr[i][0]**2) / 10.0

omega_l = 2.0 * pi / Tl
d_omega_l = omega_l * d_Tl / Tl
omega_r = 2.0 * pi / Tr
d_omega_r = omega_r * d_Tr / Tr

omega = 0.5 * (omega_l + omega_r)
d_omega = 0.5 * sqrt(d_omega_l**2 + d_omega_r**2)

omega_spec = 2.0 * pi * f
d_omega_spec = 2.0 * pi * d_f


## Output
if __name__ == '__main__':
  print(dpr.tbl([
    dpr.lst(gen.l, gen.d_l, name='l', unit='m'),
    dpr.lst(Tl, d_Tl, name='TL', unit='s'),
    dpr.lst(Tr, d_Tr, name='TR', unit='s'),
    dpr.lst(omega_l, d_omega_l, name='ω_L', unit='s'),
    dpr.lst(omega_r, d_omega_r, name='ω_R', unit='s')
  ]))

  print(dpr.tbl([
    dpr.lst(gen.l, gen.d_l, name='l', unit='m'),
    dpr.lst(omega, d_omega, name='ω', prefix=False, unit='1/s'),
    dpr.lst(omega_spec, d_omega_spec, name='ω_spec', prefix=False, unit='1/s'),
    dpr.dev(omega, d_omega, omega_spec, d_omega_spec, name='ω, ω_spec'),
  ], name='Antisymmetric oscillation frequencys'))
