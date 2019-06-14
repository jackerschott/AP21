import numpy as np
import scipy.constants as cs
from numpy import pi, sqrt

import data.print as dpr
import general as gen
import nocoupling as nc


## Data
tl = np.array([[1.78, 17.92], [2.91, 19.08], [1.84, 17.92]])
tr = np.array([[1.79, 17.93], [1.34, 17.43], [1.85, 17.88]])
d_tl = np.array([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
d_tr = np.array([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
f = np.array([0.620, 0.620, 0.620])
d_f = np.array([0.008, 0.005, 0.008])


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
omega_nc_array = np.full_like(omega, nc.omega)
d_omega_nc_array = np.full_like(omega, nc.d_omega)
if __name__ == "__main__":
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
    dpr.dev(omega, d_omega, omega_nc_array, d_omega_nc_array, name='ω, ω_nc'),
    dpr.dev(omega_spec, d_omega_spec, omega_nc_array, d_omega_nc_array, name='ω_spec, ω_nc')
  ], name='Symmetric oscillation frequencys'))
