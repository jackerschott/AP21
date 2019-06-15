import numpy as np
import scipy.constants as cs
from numpy import pi, sqrt

import dataproc.print as dpr
import general as gen
import symmetric as sym
import asymmetric as asym


## Data
tl = np.array([[0.65, 14.06], [2.75, 17.29], [2.09, 17.86]]) # 5 Periods and 10 with beats
tr = np.array([[1.10, 15.01], [2.34, 16.90], [1.71, 18.07]])
tl_b = np.array([[8.78, 52.24], [8.81, 84.62], [17.11, 187.32]])
tr_b = np.array([[4.68, 47.50], [16.08, 91.78], [33.98, 204.40]])
d_tl = np.array([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
d_tr = np.array([[0.1, 0.1], [0.1, 0.1], [0.1, 0.1]])
d_tl_b = np.array([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5]])
d_tr_b = np.array([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5]])
f = np.array([0.735, 0.686, 0.649])
d_f = np.array([0.009, 0.006, 0.003])
f_b = np.array([0.620, 0.620, 0.619])
d_f_b = np.array([0.009, 0.006, 0.003])


## Data preparation
Tl, d_Tl = np.zeros(3), np.zeros(3)
Tr, d_Tr = np.zeros(3), np.zeros(3)
Tl_b, d_Tl_b = np.zeros(3), np.zeros(3)
Tr_b, d_Tr_b = np.zeros(3), np.zeros(3)
for i in range(3):
  Tl[i] = (tl[i][1] - tl[i][0]) / 10.0
  d_Tl[i] = sqrt(d_tl[i][0]**2 + d_tl[i][1]**2) / 10.0
  Tr[i] = (tr[i][1] - tr[i][0]) / 10.0
  d_Tr[i] = sqrt(d_tr[i][0]**2 + d_tr[i][1]**2) / 10.0
  Tl_b[i] = 2.0 * (tl_b[i][1] - tl_b[i][0]) / 5.0
  d_Tl_b[i] = 2.0 * sqrt(d_tl_b[i][0]**2 + d_tl_b[i][1]**2) / 5.0
  Tr_b[i] = 2.0 * (tr_b[i][1] - tr_b[i][0]) / 5.0
  d_Tr_b[i] = 2.0 * sqrt(d_tr_b[i][0]**2 + d_tr_b[i][1]**2) / 5.0

omega_l = 2.0 * pi / Tl
d_omega_l = omega_l * d_Tl / Tl
omega_r = 2.0 * pi / Tr
d_omega_r = omega_r * d_Tr / Tr
omega_l_b = 2.0 * pi / Tl_b
d_omega_l_b = omega_l_b * d_Tl_b / Tl_b
omega_r_b = 2.0 * pi / Tr_b
d_omega_r_b = omega_r_b * d_Tr_b / Tr_b

omega = 0.5 * (omega_l + omega_r)
d_omega = 0.5 * sqrt(d_omega_l**2 + d_omega_r**2)
omega_b = 0.5 * (omega_l_b + omega_r_b)
d_omega_b = 0.5 * sqrt(d_omega_l_b**2 + d_omega_r_b**2)

omega_spec_2 = 2.0 * pi * f
d_omega_spec_2 = 2.0 * pi * d_f
omega_spec_1 = 2.0 * pi * f_b
d_omega_spec_1 = 2.0 * pi * d_f_b


## Evaluation
omega_spec = 0.5 * (omega_spec_1 + omega_spec_2)
d_omega_spec = 0.5 * sqrt(d_omega_spec_1**2 + d_omega_spec_2**2)
omega_b_spec = 0.5 * (omega_spec_2 - omega_spec_1)
d_omega_b_spec = 0.5 * sqrt(d_omega_spec_1**2 + d_omega_spec_2**2)

omega_theo = 0.5 * (sym.omega + asym.omega)
d_omega_theo = 0.5 * sqrt(sym.omega**2 + asym.d_omega**2)
omega_theo_b = 0.5 * (asym.omega - sym.omega)
d_omega_theo_b = 0.5 * sqrt(sym.d_omega**2 + asym.d_omega**2)

kappa = (asym.omega**2 - sym.omega**2) / (sym.omega**2 + asym.omega**2)
d_kappa = kappa * sqrt(4 * (asym.omega**2 * asym.d_omega**2 + sym.omega**2 * sym.d_omega**2)
                        * (1 / (asym.omega**2 - sym.omega**2)**2 + 1 / (sym.omega**2 + asym.omega**2)**2))

kappa_ratio = np.zeros(len(kappa) - 1)
d_kappa_ratio = np.zeros(len(d_kappa) - 1)
l2_ratio = np.zeros(len(gen.l) - 1)
d_l2_ratio = np.zeros(len(gen.d_l) - 1)
for i in range(len(kappa_ratio)):
  kappa_ratio[i] = kappa[i + 1] / kappa[i]
  d_kappa_ratio[i] = kappa_ratio[i] * sqrt((d_kappa[i + 1] / kappa[i + 1])**2 + (d_kappa[i] / kappa[i])**2)
  l2_ratio[i] = gen.l[i + 1]**2 / gen.l[i]**2
  d_l2_ratio[i] = l2_ratio[i] * sqrt((2 * gen.d_l[i + 1] / gen.l[i + 1])**2 + (2 * gen.d_l[i] / gen.l[i])**2)


## Output
if __name__ == '__main__':
  print(dpr.tbl([
    dpr.lst(gen.l, gen.d_l, name='l', unit='m'),
    dpr.lst(Tl_b, d_Tl_b, name='TL', unit='s'),
    dpr.lst(Tr_b, d_Tr_b, name='TR', unit='s'),
    dpr.lst(omega_l_b, d_omega_l_b, name='ω_L', unit='s'),
    dpr.lst(omega_r_b, d_omega_r_b, name='ω_R', unit='s')
  ]))

  print(dpr.tbl([
    dpr.lst(gen.l, gen.d_l, name='l', unit='m'),
    dpr.lst(Tl_b, d_Tl_b, name='TL', unit='s'),
    dpr.lst(Tr_b, d_Tr_b, name='TR', unit='s'),
    dpr.lst(omega_l_b, d_omega_l_b, name='ω_L', unit='s'),
    dpr.lst(omega_r_b, d_omega_r_b, name='ω_R', unit='s')
  ]))
  
  print(dpr.tbl([
    dpr.lst(omega, d_omega, name='ω1', prefix=False, unit='1/s'),
    dpr.lst(omega_spec, d_omega_spec, name='ω2', prefix=False, unit='1/s'),
    dpr.lst(omega_theo, d_omega_theo, name='ω_E', prefix=False, unit='1/s'),
    dpr.dev(omega, d_omega, omega_spec, d_omega_spec, name='ω1, ω2'),
    dpr.dev(omega, d_omega, omega_theo, d_omega_theo, name='ω1, ω_E'),
    dpr.dev(omega_spec, d_omega_spec, omega_theo, d_omega_theo, name='ω2, ω_E'),
  ], name='Mixed excitation frequencys'))
  print(dpr.tbl([
    dpr.lst(omega_b, d_omega_b, name='ω1', prefix=False, unit='1/s'),
    dpr.lst(omega_b_spec, d_omega_b_spec, name='ω2', prefix=False, unit='1/s'),
    dpr.lst(omega_theo_b, d_omega_theo_b, name='ω_E', prefix=False, unit='1/s'),
    dpr.dev(omega_b, d_omega_b, omega_b_spec, d_omega_b_spec, name='ω1, ω2'),
      dpr.dev(omega_b, d_omega_b, omega_theo_b, d_omega_theo_b, name='ω1, ω_E'),
    dpr.dev(omega_b_spec, d_omega_b_spec, omega_theo_b, d_omega_theo_b, name='ω2, ω_E'),
  ], name='Mixed excitation beat frequencys'))
  print()

  print(dpr.tbl([
    dpr.lst(kappa, d_kappa, name='κ', prefix=False)
  ], name='coupling factors'))

  print(dpr.tbl([
    dpr.lst(kappa_ratio, d_kappa_ratio, name='κ_ratio', prefix=False),
    dpr.lst(l2_ratio, d_l2_ratio, name='l²_ratio', prefix=False),
    dpr.dev(kappa_ratio, d_kappa_ratio, l2_ratio, d_l2_ratio, name='κ_ratio, l²_ratio')
  ], name='coupling factor ratios'))
