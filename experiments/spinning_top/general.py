from numpy import sqrt
import scipy.constants as cs

import datproc.print as dpr

import nutation_I_x as n1
import nutation_I_x_2 as n2
from precession_I_z import I_z, d_I_z

print(dpr.val(I_z / (cs.gram * cs.centi**2), d_I_z / (cs.gram * cs.centi**2),
  name='I_z', unit='g cm^2'))
print()
print(dpr.val(n1.I_x / (cs.gram * cs.centi**2), n1.d_I_x / (cs.gram * cs.centi**2),
  name='I_x', unit='g cm^2'))
print(dpr.val(n2.I_x / (cs.gram * cs.centi**2), n2.d_I_x / (cs.gram * cs.centi**2),
  name='I_x', unit='g cm^2'))
print(dpr.dev(n1.I_x, n1.d_I_x, n2.I_x, n2.d_I_x, name='I_x, I_x'))
print(dpr.dev(n1.I_x, n1.d_I_x, I_z, d_I_z, name='I_x, I_z'))
print(dpr.dev(n2.I_x, n2.d_I_x, I_z, d_I_z, name='I_x, I_z'))
