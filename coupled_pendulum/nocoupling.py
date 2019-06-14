import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cs
from numpy import pi, sqrt

import data.print as dpr

## Data
tl = np.array([1.03, 17.23])
d_tl = np.array([0.1, 0.1])

tr = np.array([1.41, 17.56])
d_tr = np.array([0.1, 0.1])


## Data preparation
Tl = (tl[1] - tl[0]) / 10.0
d_Tl = sqrt(d_tl[0]**2 + d_tl[1]**2) / 10.0

Tr = (tr[1] - tr[0]) / 10.0
d_Tr = sqrt(d_tr[0]**2 + d_tr[1]**2) / 10.0

omega_l = 2.0 * pi / Tl
d_omega_l = omega_l * d_Tl / Tl
omega_r = 2.0 * pi / Tr
d_omega_r = omega_r * d_Tr / Tr

omega = 0.5 * (omega_l + omega_r)
d_omega = 0.5 * sqrt(d_omega_l**2 + d_omega_r**2)


## Output
if __name__ == "__main__":
  print(dpr.val('T_L', Tl, d_Tl, unit='s'))
  print(dpr.val('T_R', Tr, d_Tr, unit='s'))
  print(dpr.val('ω_L', omega_l, d_omega_l, unit='1/s'))
  print(dpr.val('ω_R', omega_r, d_omega_r, unit='1/s'))
  print(dpr.val('ω', omega, d_omega, unit='1/s'))
