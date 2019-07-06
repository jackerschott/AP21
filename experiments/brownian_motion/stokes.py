import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt, exp
import os
import scipy.constants as cs
from scipy.optimize import curve_fit
from scipy.stats import norm

import datproc.plot as dpl
import datproc.print as dpr

## General
output = __name__ == '__main__'

## Data
a = 755.0 * cs.nano
d_a = 30.0 * cs.nano

T = 26.1 + cs.zero_Celsius
d_T = 0.1

eta = 7.0e-4
d_eta = 0.5e-4

t, x, y = np.loadtxt('data/brownian_motion/particle_positions.dat', usecols=(1,2,3), unpack=True)

## Data processing
x *= cs.micro
y *= cs.micro

## Evaluation
if output:
  plt.subplots(num=1)
  plt.xlabel(r'$x$ / $\mu$m')
  plt.ylabel(r'$y$ / $\mu$m')

  points = dpl.to_units(x, y, x_unit=cs.micro, y_unit=cs.micro)
  plt.plot(*points, marker='o')

delta_t = t[1:] - t[:-1]
delta_x = x[1:] - x[:-1]
delta_y = y[1:] - y[:-1]

r2 = delta_x**2 + delta_y**2

d_delta_t = np.std(delta_t) / sqrt(len(delta_t))
delta_t = np.mean(delta_t)

r2_ = np.delete(r2, 66)

d_r2_mean = np.std(r2_) / sqrt(len(r2_))
r2_mean = np.mean(r2_)

if output:
  print(dpr.val(delta_t, d_delta_t, name='Δt', unit='s'))
  print(dpr.val(r2_mean / cs.micro**2, d_r2_mean / cs.micro**2,
    name='r2', unit='μm²', prefix=False))
  print()

D = r2_mean / (4 * delta_t)
d_D = D * sqrt((d_r2_mean / r2_mean)**2 + (d_delta_t / delta_t)**2)

k = 6 * pi * eta * a * r2_mean / (4 * T * delta_t)
d_k = k * sqrt((d_eta / eta)**2 + (d_a / a)**2 + (d_r2_mean / r2_mean)**2 + (d_T / T)**2 + (d_delta_t / delta_t)**2)

if output:
  print(dpr.val(D / cs.micro**2, d_D / cs.micro**2, name='D', unit='μm^2 / s', prefix=False, exp_to_fix=0))
  print(dpr.val(k, d_k, name='k', unit='J / K', prefix=False))
  print(dpr.dev(k, d_k, cs.k, name='k, k_lit'))

## Save plots
if output:
  fig_folder_path = 'figures/brownian_motion'
  if not os.path.exists(fig_folder_path):
    os.makedirs(fig_folder_path)

  for i in plt.get_fignums():
    plt.figure(i).savefig(os.path.join(fig_folder_path, 'fig%i.pgf' % i), bbox_inches='tight', pad_inches=0.0)
    plt.figure(i).savefig(os.path.join(fig_folder_path, 'fig%i.pdf' % i), bbox_inches='tight')
