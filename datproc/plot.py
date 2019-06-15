import matplotlib.pyplot as plt
import numpy as np

from matplotlib.collections import LineCollection
from numpy import sqrt, exp
from scipy.optimize import curve_fit

# Settings
linreg_change = 0.00001 # min relative change per step to end linear regression
minfloat = 1e-80 # replaces zeros in linreg
linspace_res = 2000 # resolution for nplinspace


def initplot(num=0, nrows=1, ncols=1, title='', xlabel='', ylabel='', scale='linlin', fignum=False):
  fig, axs = plt.subplots(nrows=nrows, ncols=ncols)
  if fignum:
    if nrows != 1 or ncols != 1:
      st = fig.suptitle('Diagramm ' + str(num) + ': ' + title, fontsize='14')
    else:
      plt.title('Diagramm ' + str(num) + ': ' + title, fontsize='14')
  else:
    if nrows != 1 or ncols != 1:
      st = fig.suptitle(title, fontsize='14')
    else:
      plt.title(title, fontsize='14')
  for ax in np.array([axs]).reshape(-1):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, which='both')
    if (scale == 'linlin'):
      ax.ticklabel_format(style='sci', axis='both', scilimits=(-2,3))
    elif (scale == 'linlog'):
      ax.set_yscale('log')
      ax.ticklabel_format(style='sci', axis='x', scilimits=(-2,3))
    elif (scale == 'loglin'):
      ax.set_xscale('log')
      ax.ticklabel_format(style='sci', axis='y', scilimits=(-2,3))
    elif (scale == 'loglog'):
      ax.set_yscale('log')
      ax.set_xscale('log')
  fig.set_size_inches(11.69,8.27)
  fig.tight_layout()
  if nrows != 1 or ncols != 1:
    st.set_y(0.97)
    fig.subplots_adjust(top=0.90)

def set_axis(num):
  plt.sca(plt.gcf().axes[num])

def plot(x, y, label='', color=None):
  plot = plt.plot(x, y, label=label, color=color)
  if label != '':
    plt.legend()
  return plot

def plotdata(x, y, dy=None, dx=None, label='', color=None, connect=False):
  plot = plt.errorbar(x=x, y=y, yerr=dy, xerr=dx, label=label, color=color, fmt='o', markersize=3, capsize=5)
  # Other plot options
  for cap in plot[1]:
    cap.set_markeredgewidth(1)
  if (connect == True):
    if (color == None):
      color = plot[0].get_color()
    plt.plot(x, y, color=color)
  if label != '':
    plt.legend()
  return plot

def linreg(x, y, dy, dx=None, fit_range=None, plot=False, graphname='', legend=False, scaleReg=False):
  # Set fit range if None
  if (fit_range == None):
    fit_range = range(len(x))
  
  # Regression iteration, for dx is None only one iteration is needed
  def linreg_iter(x, y, dy):
    [s0, s1, s2, s3, s4] = [0.0, 0.0, 0.0, 0.0, 0.0]
    for i in fit_range:
      if (dy[i] == 0.0):
        dy[i] = minfloat
      s0 += dy[i]**-2
      s1 += x[i] * dy[i]**-2
      s2 += y[i] * dy[i]**-2
      s3 += x[i]**2 * dy[i]**-2
      s4 += x[i] * y[i] * dy[i]**-2
    eta = s0 * s3 - s1**2
    g = (s0 * s4 - s1 * s2) / eta
    dg = np.sqrt(s0 / eta)
    b = (s3 * s2 - s1 * s4) / eta
    db = np.sqrt(s3 / eta)
    return [g, dg, b, db]

  # Compute slope and axis intercept
  iter0 = linreg_iter(x, y, dy)
  result = []
  dx_ = dx
  if (dx is None):
    dx_ = np.zeros(len(x))
    result = iter0
  else:
    g = iter0[0]
    g_old = g * (1 - 2 * linreg_change)
    dy_ = dy
    while (abs(1 - g_old / g) >= linreg_change):
      g_old = g
      dy_ = np.sqrt((g * dx)**2 + dy_**2)
      g = linreg_iter(x, y, dy_)[0]
    result = linreg_iter(x, y, dy_)
  
  # Plot
  if (plot):
    # Get data for regression line plot
    [g, dg, b, db] = result
    min_x = np.argmin(x)
    max_x = np.argmax(x)
    xint = np.linspace(x[min_x] - dx_[min_x], x[max_x] + dx_[max_x], linspace_res)
    yfit = g * xint + b
    yerr = (g + dg) * xint + (b - db)

    # Plot data points
    data_plot = plotdata(x=x, y=y, dy=dy, dx=dx, label=graphname)

    # Plot regression line and uncertainty line
    color = data_plot[0].get_color()
    ax = plt.gca()
    if scaleReg:
      plt.plot(xint, yfit, marker='', color=color)
      plt.plot(xint, yerr, marker='', linestyle='dashed', color=color)
    else:
      reg_line = LineCollection([np.column_stack((xint, yfit))], colors=color)
      reg_err_line = LineCollection([np.column_stack((xint, yerr))], colors=color, linestyles='dashed')
      ax.add_collection(reg_line, autolim=False)
      ax.add_collection(reg_err_line, autolim=False)

    # Add legend
    if (legend):
      plt.legend(['Fit', 'Fit uncertainty'])
    elif (graphname != ''):
      plt.legend()
  return result

def expreg(x, y, dy, dx=None, fit_range=None, plot=True):
  if fit_range == None:
    fit_range = range(len(x))
  expo, dexpo, _yitc, _dyitc = linreg(x, np.log(y), dy/y, dx, fit_range=fit_range)
  yitc = exp(_yitc)
  dyitc = yitc * _dyitc
  result = [expo,dexpo,yitc,dyitc]
  if (plot):
    dx_ = dx
    if dx == None:
      dx_ = np.zeros(len(x))
    min_x = np.argmin(x)
    max_x = np.argmax(x)
    xint = np.linspace(x[min_x] - dx_[min_x], x[max_x] + dx_[max_x], 1000)
    yfit = yitc * exp(expo * xint)
    yerr = (yitc - dyitc) * exp((expo + dexpo) * xint)
    data_plot = plotdata(x=x, y=y, dy=dy, dx=dx)
    color = data_plot[0].get_color()
    left, right = plt.xlim()
    top, bottom = plt.ylim()
    plt.plot(xint, yfit, marker='', color=color)
    plt.plot(xint, yerr, marker='', linestyle='dashed', color=color)
    plt.xlim(left, right)
    plt.ylim(top, bottom)
  return result

def fit(x, y, dy, f, p0=None, fit_range=None, plot=True):
  if fit_range == None:
    fit_range = range(len(x))
  x_fit = [x[i] for i in fit_range]
  y_fit = [y[i] for i in fit_range]
  dy_fit = [dy[i] for i in fit_range]
  p, d_p = curve_fit(f, x_fit, y_fit, sigma=dy_fit, p0=p0)
  if plot:
    xint = np.linspace(np.min(x), np.max(x), linspace_res)
    yfit = f(xint, *p)
    data_plot = plotdata(x, y, dy)
    color = data_plot[0].get_color()
    left, right = plt.xlim()
    top, bottom = plt.ylim()
    plt.plot(xint, yfit, marker='', color=color)
    plt.xlim(left, right)
    plt.ylim(top, bottom)
  return (p, sqrt(np.diag(d_p)))

def lin_yerr(x, dx, y, dy):
    g = linreg(x, y, dx, dy)
    new_dy = [np.sqrt(dy[i]**2 + (g * dx[i])**2) for i in range(len(dy))]
    return new_dy


def savefigs(path):
  # Save figures in 'path' as figN.pdf, where N is the figures number
  for i in plt.get_fignums():
    plt.figure(i).savefig(path + '/fig' + str(i) +'.pdf', papertype='a4', orientation='landscape', bbox_inches='tight', pad_inches=0.6, format='pdf')
