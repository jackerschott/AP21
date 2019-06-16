import datproc.plot as dp
import datproc.print as dpr
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cs
from numpy import pi, sqrt, sin, cos, tan, sinc
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.signal import find_peaks



## General

print()

fa = lambda a: np.array(a, dtype=float)

# The cut is performed after cut as index. lshift = 0 returns original array.
def dat_overlap(arr, d_arr, cut, lshift):
  c = cut
  s = lshift
  arr_ = np.concatenate((arr[:c-s], np.mean([arr[c-s:c], arr[c:c+s]], axis=0)))
  arr_ = np.concatenate((arr_, arr[c+s:]))
  d_arr_ = np.concatenate((d_arr[:c-s], sqrt(d_arr[c-s:c]**2 + d_arr[c:c+s]**2) / 2))
  d_arr_ = np.concatenate((d_arr_, d_arr[c+s:]))
  return arr_, d_arr_

plotTitles = [
  r'Lineare Regression der Analysierspaltbreite $d$ in Abhängigkeit der Position $x$' '\n'
    r'der jeweiligen gerade so verdeckten Minima.',
  r'Positionen $x$ der Minima und Maxima des Einzelspaltes in Abhängigkeit der Ordnung $n$.',
  r'Intensität $I$ des Doppelspalt-Fourierbilds in Abhängigkeit der Beugungsordnung $n$' '\n'
    r'mit dem Fourierbild des Einzelspalts als Referenz',
  r'Intensität $I$ des Einzelspalt-Objektbilds in Abhängigkeit der Position $y$' '\n'
    r'unter Beachtung von einer bis drei Beugungsordnungen',
  r'Intensität $I$ des Doppelspalt-Objektbilds in Abhängigkeit der Position $y$' '\n'
    r'im Fall (a) (doppeltes Gaußprofil) und im Fall (b) (Plateau)'
]



## Measured data

# General
lda = 635 * cs.nano
f = 80 * cs.milli
px = 14 * cs.micro


# Maxima, minima of the single slit fourier image
n_sf = fa([0, 1, 1, 2, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6])
I_ug_sf = 72.0
d_I_ug_sf = 14.00
I_sf = fa([3700, 90, 290, 84, 3340, 130, 1330, 110, 640, 100, 410, 96, 287, 99])
d_I_sf = fa([50, 10, 20, 10, 50, 20, 40, 20, 30, 20, 30, 10, 30, 10])
x_sf = fa([1073, 992, 955, 913, 957, 912, 877, 834, 794, 755, 717, 676, 639, 598])
d_x_sf = fa([5] * len(x_sf))


# Maxima, minima of the double slit fourier image
n_df = fa([0, 1, 1, 2, 1, 2, 2, 3, 3, 4])
I_ug_df = 72.0
d_I_ug_df = 14.0
I_df = fa([3530, 120, 2230, 90, 3590, 100, 390, 90, 435, 90])
d_I_df = fa([50, 30, 50, 20, 50, 20, 40, 20, 30, 20])
x_df = fa([1219, 1202, 1188, 1169, 1184, 1164, 1155, 1140, 1114, 1098])
d_x_df = fa([5, 5, 5, 5, 5, 5, 5, 30, 5, 5])


# Single slit object image
I_so = np.array([fa([1580]), fa([1080, 1580]), fa([1440, 1120, 1530])])
I_ug_so = fa([90, 86, 89])
d_I_so = np.array([fa([40]), fa([20, 20]), fa([30, 30, 50])])
d_I_ug_so = fa([3, 5, 5])

x1_so = 939
d_x1_so = 3
x2_so = 1070
d_x2_so = 3
b_so = 682 * cs.milli
d_b_so = 5 * cs.milli


# Double slit object image
x1_do = 714
d_x1_do = 15
x2_do = 819
d_x2_do = 15
x3_do = 1000
d_x3_do = 15
x4_do = 1105
d_x4_do = 15
b_do = 137 * cs.milli
d_b_do = 3 * cs.milli

w_a_do = 0.23 * cs.milli
d_w_a_do = 0.01 * cs.milli
w_b_do = 0.12 * cs.milli
d_w_b_do = 0.01 * cs.milli


# Abscissa calibration (single slit)
n_sc = fa([1, 2, 3, 4])
d_sc = 2 * fa([0.235, 0.440, 0.650, 0.860]) * cs.milli
d_d_sc = 2 * fa([0.005, 0.005, 0.005, 0.005]) * cs.milli



## Data Preparation


# Maxima, minima of the single slit fourier image
I_sf -= I_ug_sf
d_I_sf = sqrt(d_I_sf**2 + d_I_ug_sf**2)
I_sf_ = np.concatenate((
  I_sf[:4] / I_sf[0],
  I_sf[4:] / I_sf[4] * I_sf[2] / I_sf[0]
))
d_I_sf_ = np.concatenate((
  I_sf_[:4] * sqrt((d_I_sf[:4] / I_sf[:4])**2 + (d_I_sf[0] / I_sf[0])**2),
  I_sf_[4:] * sqrt((d_I_sf[4:] / I_sf[4])**2 + (d_I_sf[4] / I_sf[4])**2 + (d_I_sf[2] / I_sf[2])**2 + (d_I_sf[0] / I_sf[0])**2)
))

n_sf = np.concatenate((n_sf[:4], n_sf[6:]))
I_sf, d_I_sf = dat_overlap(I_sf_, d_I_sf_, 4, 2)
x_sf, d_x_sf = dat_overlap(x_sf, d_x_sf, 4, 2)

n_min_sf = n_sf[1::2]
I_min_sf = I_sf[1::2]
d_I_min_sf = d_I_sf[1::2]
x_min_sf = x_sf[1::2]
d_x_min_sf = d_x_sf[1::2]

I_max_sf = I_sf[::2]
d_I_max_sf = d_I_sf[::2]
x_max_sf = x_sf[::2]
d_x_max_sf = d_x_sf[::2]


# Maxima, minima of the double slit fourier image
I_df -= I_ug_df
d_I_df = sqrt(d_I_df**2 + d_I_ug_df**2)
I_df_ = np.concatenate((
  I_df[:4] / I_df[0],
  I_df[4:] / I_df[4] * I_df[2] / I_df[0]
))
d_I_df_ = np.concatenate((
  I_df_[:4] * sqrt((d_I_df[:4] / I_df[:4])**2 + (d_I_df[0] / I_df[0])**2),
  I_df_[4:] * sqrt((d_I_df[4:] / I_df[4:])**2 + (d_I_df[4] / I_df[4])**2 + (d_I_df[2] / I_df[2])**2 + (d_I_df[0] / I_df[0])**2)
))

n_df = np.concatenate((n_df[:4], n_df[6:]))
I_df, d_I_df = dat_overlap(I_df_, d_I_df_, 4, 2)
x_df, d_x_df = dat_overlap(x_df, d_x_df, 4, 2)

n_min_df = n_df[1::2]
I_min_df = I_df[1::2]
d_I_min_df = d_I_df[1::2]
x_min_df = x_df[1::2]
d_x_min_df = d_x_df[1::2]

I_max_df = I_df[::2]
d_I_max_df = d_I_df[::2]
x_max_df = x_df[::2]
d_x_max_df = d_x_df[::2]


# Single slit object image
for i in range(len(I_so)):
  I_i_max_so = np.argmax(I_so[i])
  d_I_so[i] = I_so[i] / I_so[i][I_i_max_so] * sqrt((d_I_so[i] / I_so[i])**2 + (d_I_so[i][I_i_max_so] / I_so[i][I_i_max_so])**2)
  I_so[i] = I_so[i] / I_so[i][I_i_max_so]

wpx_so = x2_so - x1_so
d_wpx_so = sqrt(d_x1_so**2 + d_x2_so**2)
M_so = b_so / f - 1
d_M_so = d_b_so / f
w_so = wpx_so * px / M_so
d_w_so = w_so * sqrt((d_wpx_so / wpx_so)**2 + (d_M_so / M_so)**2)

print(dpr.val('M', M_so, d_M_so))
print(dpr.val('w_px', wpx_so, d_wpx_so, prefix=False, unit='px'))
print(dpr.val('w', w_so, d_w_so, unit='m'))
print()


# Double slit object image
wpx1_do = x2_do - x1_do
d_wpx1_do = sqrt(d_x1_do**2 + d_x2_do**2)
wpx2_do = x4_do - x3_do
d_wpx2_do = sqrt(d_x3_do**2 + d_x4_do**2)

gpx_do = x3_do - x2_do + wpx1_do / 2 + wpx2_do / 2
d_gpx_do = sqrt(d_x2_do**2 + d_x3_do**2 + (d_wpx1_do**2 + d_wpx2_do**2) / 4)

M_do = b_do / f - 1
d_M_do = d_b_do / f
w1_do = wpx1_do * px / M_do
d_w1_do = w1_do * sqrt((d_wpx1_do / wpx1_do)**2 + (d_M_do / M_do)**2)
w2_do = wpx2_do * px / M_do
d_w2_do = w2_do * sqrt((d_wpx2_do / wpx2_do)**2 + (d_M_do / M_do)**2)
w_mean_do = (w1_do + w2_do) / 2
d_w_mean_do = sqrt(d_w1_do**2 + d_w2_do**2) / 2

g_do = gpx_do * px / M_do
d_g_do = g_do * sqrt((d_gpx_do / gpx_do)**2 + (d_M_do / M_do)**2)

print(dpr.val('M', M_do, d_M_do))
print(dpr.val('w1_px', wpx1_do, d_wpx1_do, prefix=False, unit='px'))
print(dpr.val('w2_px', wpx2_do, d_wpx2_do, prefix=False, unit='px'))
print(dpr.val('g_px', gpx_do, d_gpx_do, prefix=False, unit='px'))
print(dpr.val('w1', w1_do, d_w1_do, unit='m'))
print(dpr.val('w2', w2_do, d_w2_do, unit='m'))
print(dpr.val('w', w_mean_do, d_w_mean_do, unit='m'))
print(dpr.val('g', g_do, d_g_do, unit='m'))
print()



## Evaluation

# Fourier and object image functions
def I_slit(n):
  return sinc(n)**2
def d_I_slit(n, d_n):
  return 2 * abs((sinc(n) - cos(pi * n)) * sinc(n) * d_n / n)
def I_dslit(n, v):
  return (sinc(n) * cos(pi * v * n))**2
def slit_mod_kernel(n, y):
  return 2 * sinc(n) * cos(2 * pi * n * y)
def dslit_mod_kernel(n, y, g):
  return 4 * sinc(n) * cos(pi * n * g) * cos(2 * pi * n * y)


# Abscissa calibration (single slit)
i_min_sc = np.array(n_sc - 1, dtype=int)
dp.initplot(num=2, title=plotTitles[0], xlabel=r'$x$ / px', ylabel=r'$d$ / mm', fignum=True)
s_sc, d_s_sc, b_sc, d_b_sc = dp.linreg(x_min_sf[i_min_sc], d_sc / cs.milli, d_d_sc / cs.milli, d_x_min_sf[i_min_sc], plot=True)
s_sc *= cs.milli
d_s_sc *= cs.milli
b_sc *= cs.milli
d_b_sc *= cs.milli

print(dpr.val('s', s_sc, d_s_sc, unit='m / px'))
print(dpr.val('b', b_sc, d_b_sc, unit='m'))
print()


# 1. Single slit fourier image
dp.initplot(num=1, title=plotTitles[1], xlabel=r'$n$', ylabel=r'$x$ / px', fignum=True)
s_sf, d_s_sf, b_sf, d_b_sf = dp.linreg(n_min_sf, x_min_sf, d_x_min_sf, plot=True)
n_max_sf = (x_max_sf - b_sf) / s_sf
d_n_max_sf = abs(n_max_sf) * sqrt((d_x_max_sf**2 + d_b_sf**2) / (x_max_sf - b_sf)**2 + (d_s_sf / s_sf)**2)
dp.plotdata(n_max_sf, x_max_sf, d_x_max_sf, d_n_max_sf)

w_sf = 2 * lda * f / (s_sc * s_sf)
d_w_sf = w_sf * sqrt((d_s_sc / s_sc)**2 + (d_s_sf / s_sf)**2)

def tan_id(x):
  return tan(x) - x

n_max_sf_theo = np.zeros_like(n_max_sf)
for i in range(len(n_max_sf_theo)):
  n_max_sf_theo[i] = brentq(tan_id, -pi/2 + i * pi + pi / 10000, pi/2 + i * pi - pi / 10000) / pi

I_max_sf_theo = I_slit(n_max_sf)
d_I_max_sf_theo = d_I_slit(n_max_sf, d_n_max_sf)

print(dpr.val('s', s_sf, d_s_sf, prefix=False, unit='px'))
print(dpr.val('b', b_sf, d_b_sf, prefix=False, unit='px'))
print()
print(dpr.tbl([
  dpr.lst(n_max_sf, d_n_max_sf, name='n_o', expToFix=0),
  dpr.lst(n_max_sf_theo, name='n_t', expToFix=0),
  dpr.sig('n_o, n_t', n_max_sf, d_n_max_sf, n_max_sf_theo, perc=True)
]))
print(dpr.tbl([
  dpr.lst(I_max_sf, d_I_max_sf, name='I_o', unit='I0', prefix=False, expToFix=0),
  dpr.lst(I_max_sf_theo, d_I_max_sf_theo, name='I_t', unit='I0', prefix=False, expToFix=0),
  dpr.sig('I_o, I_t', I_max_sf, d_I_max_sf, I_max_sf_theo, d_I_max_sf_theo, perc=True)
]))
print()


# 2. Double slit fourier image
n_df = np.linspace(-2, 2, 1000)
v_df = g_do / w_mean_do
d_v_df = v_df * sqrt((d_g_do / g_do)**2 + (d_w_mean_do / w_mean_do)**2)

n_max_df_theo = fa([0, 1, 2, 4]) / v_df
d_n_max_df_theo = n_max_df_theo * d_v_df / v_df
I_max_df_theo = I_dslit(n_max_df_theo, v_df)

dp.initplot(num=3, title=plotTitles[2], xlabel=r'$n$', ylabel=r'$I$ / b.E.', fignum=True)
dp.plot(n_df, I_slit(n_df), label='Einzelspalt')
dp.plot(n_df, I_dslit(n_df, v_df), label='Doppelspalt')

print(dpr.tbl([
  dpr.lst(I_max_df, d_I_max_df, name='I_o', unit='I0', prefix=False, expToFix=0),
  dpr.lst(I_max_df_theo, name='I_t', unit='I0', prefix=False, expToFix=0),
  dpr.sig('I_o, I_t', I_max_df, d_I_max_df, I_max_df_theo, perc=True)
]))
print()


# 3. Single slit object image
dp.initplot(num=4, nrows=2, ncols=2, title=plotTitles[3], xlabel=r'$y$ / $d$', ylabel=r'$I$ / b.E.', fignum=True)
I_so_theo = np.zeros_like(I_so)
for i in range(len(I_so)):
  y_so = np.linspace(-1, 1, 100)

  I_cont_so = np.array([quad(lambda n: slit_mod_kernel(n, y), 0, i + 1) for y in y_so])
  I_cont_so = np.array([x[0]**2 for x in I_cont_so])
  I_cont_so = I_cont_so / np.max(I_cont_so)

  # Super ugly code to perform the super ugly task of finding the peaks
  I_i_max_so = find_peaks(I_cont_so)[0]
  I_i_min_so = find_peaks(-I_cont_so)[0]
  I_i_so = np.concatenate((I_i_max_so, I_i_min_so))
  I_peaks_so = I_cont_so[I_i_so]
  I_peaks_so = np.unique(np.round(I_peaks_so[I_peaks_so>0.01], decimals=10))
  if i == 2:
    I_peaks_so[0], I_peaks_so[1] = I_peaks_so[1], I_peaks_so[0]
  I_so_theo[i] = I_peaks_so

  dp.set_axis(i)
  dp.plot(y_so, I_cont_so)

for i in range(len(I_so)):
  print(dpr.tbl([
    dpr.lst(I_so[i], d_I_so[i], name='I_exp', unit='I_max', prefix=False, expToFix=0),
    dpr.lst(I_so_theo[i], name='I_theo', unit='I_max', prefix=False, expToFix=0),
    dpr.sig('I_exp, I_theo', I_so[i], d_I_so[i], I_so_theo[i], perc=True)
  ]))


# 4. Double slit object image
g_do = 2
n_do_a = 1.0
n_do_b = 0.335 # determined by trial and error
y_do = np.linspace(-2, 2, 1000)

I_do_a = np.array([quad(lambda n: dslit_mod_kernel(n, y, g_do), 0, n_do_a) for y in y_do])
I_do_a = np.array([x[0]**2 for x in I_do_a])
I_do_a = I_do_a / np.max(I_do_a)

I_do_b = np.array([quad(lambda n: dslit_mod_kernel(n, y, g_do), 0, n_do_b) for y in y_do])
I_do_b = np.array([x[0]**2 for x in I_do_b])
I_do_b = I_do_b / np.max(I_do_b)

dp.initplot(num=5, ncols=2, title=plotTitles[4], xlabel=r'$y$ / $d$', ylabel=r'$I$ / b.E.', fignum=True)
plt.ylim(0, 1)
dp.set_axis(0)
dp.plot(y_do, I_do_a)
dp.set_axis(1)
dp.plot(y_do, I_do_b)

n_do_exp_b = w_b_do * w_mean_do / (2 * f * lda)
d_n_do_exp_b = n_do_exp_b * sqrt((d_w_b_do / w_b_do)**2 + (d_w_mean_do / w_mean_do)**2)

k_do_b = 2 * pi / w_mean_do * n_do_b
k_do_exp_b = 2 * pi / w_mean_do * n_do_exp_b
d_k_do_exp_b = k_do_exp_b * sqrt((d_w_mean_do / w_mean_do)**2 + (d_n_do_exp_b / n_do_exp_b)**2)

print(dpr.val('n_exp', n_do_exp_b, d_n_do_exp_b))
print(dpr.val('n_theo', n_do_b))
print(dpr.sig('n_exp, n_theo', n_do_exp_b, d_n_do_exp_b, n_do_b, perc=True))
print()
print(dpr.val('k_exp', k_do_exp_b * cs.milli, d_k_do_exp_b * cs.milli, unit='1/mm'))
print(dpr.val('k_theo', k_do_b * cs.milli, unit='1/mm'))
print(dpr.sig('k_exp, k_theo', k_do_exp_b, d_k_do_exp_b, k_do_b, perc=True))
print()


# Print remaining values
print(dpr.val('w', w_sf, d_w_sf, unit='m'))
print(dpr.val('w', w_so, d_w_so, unit='m'))
print(dpr.sig('w', w_sf, d_w_sf, w_so, d_w_so, perc=True))
print()



## Show plots
dp.savefigs('figures/fourier')
