import numpy as np
from numpy import sqrt
import scipy.constants as cs

import datproc.print as dpr

# General
pL1 = 999.4 * cs.hecto
pL2 = 999.5 * cs.hecto

kappa_L = 1.40
kappa_Ar = 1.67

# Clement, Desormes
h1 = np.array([126.0, 22.0, 137.0, 30.0, 133.0]) * cs.milli
d_h1 = np.array([1.4] * 5) * cs.milli
h3 = np.array([22.0, 4.0, 30.0, 7.0, 28.0]) * cs.milli
d_h3 = np.array([1.4] * 5) * cs.milli

# Ruechart
V = 5460 * cs.centi**3
d_V = 5 * cs.centi**3
m = 26.006 * cs.gram
d_m = 0.002 * cs.gram
r = 15.97 / 2 * cs.milli
d_r = 0.05 / 2 * cs.milli

N = 50
TL1 = 23.1 + cs.zero_Celsius
TL2 = 23.3 + cs.zero_Celsius
tL = 46.42 / N # 47.42 / N
d_tL = 0.5 / N
TAr1 = 23.4 + cs.zero_Celsius
TAr2 = 23.5 + cs.zero_Celsius
tAr = 46.31 / N
d_tAr = 0.5 / N

# Evaluation
kappa_cd = h1 / (h1 - h3)
d_kappa_cd = kappa_cd * sqrt((1 / h1 - 1 / (h1 - h3))**2 * d_h1**2 + d_h3**2 / (h1 - h3)**2)

print(dpr.val(tL, d_tL, name='TL', unit='s', prefix=False, exp_to_fix=0))
print(dpr.val(tAr, d_tAr, name='TAr', unit='s', prefix=False, exp_to_fix=0))
print()

TL = 0.5 * (TL2 + TL1)
d_TL = 0.5 * (TL2 - TL1)
TAr = 0.5 * (TAr2 + TAr1)
d_TAr = 0.5 * (TAr2 - TAr1)
pL = 0.5 * (pL2 + pL1)
d_pL = 0.5 * (pL2 - pL1)

print(dpr.val(pL, d_pL, name='pL', unit='Pa'))
print()

kappa_rL = 4 * m * V / (r**4 * tL**2 * pL)
d_kappa_rL = kappa_rL * sqrt((d_m / m)**2 + (d_V / V)**2 + (4 * d_r / r)**2 + (2 * d_tL / tL)**2 + (d_pL / pL)**2)
kappa_rAr = 4 * m * V / (r**4 * tAr**2 * pL)
d_kappa_rAr = kappa_rAr * sqrt((d_m / m)**2 + (d_V / V)**2 + (4 * d_r / r)**2 + (2 * d_tAr / tAr)**2 + (d_pL / pL)**2)

print(dpr.tbl([
  dpr.lst(kappa_cd, d_kappa_cd, name='ĸ_cd'),
  dpr.dev(kappa_cd, d_kappa_cd, kappa_rL, d_kappa_rL, name='κ_cd, κ_r')
]))

kappa_cd = np.mean(kappa_cd)
d_kappa_cd = sqrt(np.sum(d_kappa_cd**2)) / len(d_kappa_cd)

print(dpr.val(kappa_cd, d_kappa_cd, name='ĸ_cd'))
print(dpr.val(kappa_rL, d_kappa_rL, name='ĸ_rL'))
print()
print(dpr.dev(kappa_cd, d_kappa_cd, kappa_L, name='ĸ_cd, ĸ_L'))
print(dpr.dev(kappa_rL, d_kappa_rL, kappa_L, name='ĸ_rL, ĸ_L'))
print(dpr.dev(kappa_cd, d_kappa_cd, kappa_rL, d_kappa_rL, name='ĸ_cd, ĸ_rL'))
print()
print(dpr.val(kappa_rAr, d_kappa_rL, name='ĸ_rAr'))
print()
print(dpr.dev(kappa_rAr, d_kappa_rAr, kappa_Ar, name='ĸ_rAr, ĸ_Ar'))