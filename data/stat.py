import numpy as np
from numpy import sqrt

def mv(x):
  s = 0.0
  for i in range(len(x)):
    s += x[i]
  return s / len(x)

def dsto(x):
  s = 0.0
  for i in range(len(x)):
    s += (x[i] - mv(x))**2
  return sqrt(s / (len(x) - 1))

def dsto_mv(x):
  return dsto(x) / sqrt(len(x))

def dsys_mv(x):
  return sqrt(np.sum(x**2)) / len(x)

def dtot(dsys, dsto):
  return sqrt(dsys**2 + dsto**2)

def chi2(yo, dyo, ye, dye=[]):
  if (dye == []):
    dye = [0.0 for i in range(len(ye))]
  chi2 = 0.0
  for i in range(len(yo)):
    chi2 += (yo[i] - ye[i])**2 / (dyo[i]**2 + dye[i]**2)
  return chi2

def chi2_red(yo, dyo, ye, dye=[], dof=0):
  if (dof == 0):
    dof = len(ye)
  return chi2(yo, dyo, ye, dye) / dof

def dev(yo, dyo, ye, dye=None):
  if dye is None:
    dye = np.zeros_like(dyo)
  return np.abs(yo - ye) / sqrt(dyo**2 + dye**2)