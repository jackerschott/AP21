import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt, exp
import os
import scipy.constants as cs

import datproc.plot as dpl
import datproc.print as dpr

from stokes import delta_x, delta_y

## General
output = __name__ == '__main__'

def normpdf(x, mu, sigma):
  return exp(-0.5 * (x - mu)**2 / sigma**2) / sqrt(2 * pi * sigma**2)

delta_s = np.append(delta_x, delta_y)
delta_s_ = np.delete(delta_s, 66)
if output:
  print(dpr.val(delta_s[66], name='Δs[66]', unit='m'))

delta_s_mu = np.mean(delta_s_)
delta_s_sigma = np.std(delta_s_)
if output:
  print(dpr.val(delta_s_mu, name='μ(Δs)', unit='m'))
  print(dpr.val(delta_s_sigma, name='σ(Δs)', unit='m'))

if output:
  x_gauss = np.linspace(np.min(delta_s_), np.max(delta_s_), 1000)
  y_gauss = normpdf(x_gauss / cs.micro, delta_s_mu / cs.micro, delta_s_sigma / cs.micro)

  plt.subplots(num=2)
  plt.xlabel(r'$\Delta s$ / $\mu$m')
  plt.ylabel(r'Relative Häufigkeit')
  plt.hist(delta_s_ / cs.micro, bins=20, density=True)
  plt.plot(x_gauss / cs.micro, y_gauss)


## Save plots
if output:
  fig_folder_path = 'figures/brownian_motion'
  if not os.path.exists(fig_folder_path):
    os.makedirs(fig_folder_path)

  for i in plt.get_fignums():
    plt.figure(i).savefig(os.path.join(fig_folder_path, 'fig%i.pgf' % i), bbox_inches='tight', pad_inches=0.0)
    plt.figure(i).savefig(os.path.join(fig_folder_path, 'fig%i.pdf' % i), bbox_inches='tight')
