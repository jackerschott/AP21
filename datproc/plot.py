import numpy as np
from os import path

def __to_float_array(x):
  if type(x) is list or type(x) is tuple:
    x = np.asarray(x)
  elif type(x) is not np.ndarray:
    raise TypeError('x is not an array_like object')
  
  if x.dtype != 'float' and x.dtype != 'int':
    raise TypeError('The values of x must be convertable to float')
  
  return x

def to_units(x, y, dy=None, dx=None, x_unit=1, y_unit=1):
  if dx is None and dy is None:
    return x / x_unit, y / y_unit
  if dx is None:
    return x / x_unit, y / y_unit, dy / y_unit
  if dy is None:
    return x / x_unit, y / y_unit, dx / x_unit
  else:
    return x / x_unit, y / y_unit, dy / y_unit, dx / x_unit

def linreg(x, y, dy, dx=None):
  change_max = 0.00001
  zero_replacement = 1e-80

  # Regression iteration, for dx is None only one iteration is needed
  def linreg_iter(x, y, dy):
    dy[dy == 0] = zero_replacement
    s0 = np.sum(1 / dy**2)
    s1 = np.sum(x / dy**2)
    s2 = np.sum(y / dy**2)
    s3 = np.sum(x**2 / dy**2)
    s4 = np.sum(x * y / dy**2)

    eta = s0 * s3 - s1**2
    s = (s0 * s4 - s1 * s2) / eta
    ds = np.sqrt(s0 / eta)
    b = (s3 * s2 - s1 * s4) / eta
    db = np.sqrt(s3 / eta)
    return s, ds, b, db

  x = __to_float_array(x)
  y = __to_float_array(y)
  dy = __to_float_array(dy)

  # Compute slope and axis intercept for dx not specified
  if dx is None:
    return linreg_iter(x, y, dy)

  dx = __to_float_array(dx)

  # Compute slope and axis intercept for dx specified
  dy_ = dy.copy()
  s, ds, b, db = linreg_iter(x, y, dy_)
  while True:
    s_old = s
    dy_ = np.sqrt((s * dx)**2 + dy_**2)
    s, *_ = linreg_iter(x, y, dy_)
    if abs(1 - s_old / s) < change_max:
      break
  return linreg_iter(x, y, dy_)

def x_fit_like(a, num=1000):
  return np.linspace(a[0], a[-1], num)

def linreg_lines(x, s, ds, b, db, align='flat'):
  if align != 'flat' and align != 'steep':
    raise ValueError('Unknown option for align')

  y_line = s * x + b

  flat_is_plus = abs(s - ds) > abs(s + ds)
  if (align == 'flat' and flat_is_plus
    or align == 'steep' and not flat_is_plus):
    y_uline = (s + ds) * x + (b - db)
  elif (align == 'flat' and not flat_is_plus
    or align == 'steep' and flat_is_plus):
    y_uline = (s - ds) * x + (b + db)
  
  return y_line, y_uline

def get_fig_paths(folder_path, fignums, format='pdf'):
  fig_paths = [''] * len(fignums)
  for i, num in enumerate(fignums):
    file_name = 'fig%i.%s' % (num, format)
    fig_paths[i] = path.join(folder_path, file_name)
  return fig_paths
