"""This stock site frequency plugin employs a double exponential as a simple site frequency model.

f = a1 * exp(-p * k1) + a2 * exp(-p * k2)   (f is normalized such that sum(f) = 1.0)
"""
import numpy as np

import logging
logger = logging.getLogger(__name__)

__example_param_text = """
{
  "double_exp": {
    "a1": 1.0,
    "k1": 20.0,
    "a2": 0.1,
    "k2": 5.0,
    "p0": 0.001,
    "p1": 0.2,
    "bin_cnt": 30
  }
}
"""

_description = """
This is a simple site frequency plugin. A typical parameter set resembles
""" + __example_param_text

_example_params = eval(__example_param_text)


class Model:
  def __init__(self, a1=1.0, k1=20.0, a2=0.1, k2=5.0, p0=0.001, p1=0.2, bin_cnt=21):
    """Simple double exponential model for site frequency spectrum

    :param a1: the amplitude of the fast component
    :param k1: the fast exponential constant
    :param a2: the amplitude of the slow constant
    :param k2: the slow exponential constant
    :param p0: the lower range of probability values
    :param p1: the upper range of probability values
    :param bin_cnt:  the number of bins for the frequency spectrum
    """
    p = np.linspace(p0, p1, bin_cnt)
    f = a1 * np.exp(-p * k1) + a2 * np.exp(-p * k2)
    f = f / f.sum()
    self.p, self.f = p, f

  def get_p_effective(self, p):
    """Given the single sample per-base mutation probability give us the effective per-base mutation probability given
    our site frequency spectrum

    :param p: sample per-base mutation probability
    :return:
    """
    return p / (self.p * self.f).sum()

  def get_spectrum(self):
    """Return the site frequency spectrum as a tuple of (p, f) """
    return self.p, self.f

  def __repr__(self):
    """ASCII plot of function!"""

    # We plot it as a sideways bar-graph
    size_x = 80  # columns
    scaling_factor = float(size_x) / self.f.max()
    rep_str = ''
    rep_str += ' p                 f\n'
    for p_i, f_i in zip(self.p, self.f):
      rep_str += '{:1.3f} '.format(p_i) + '-' * int(scaling_factor * f_i + 0.5) + '{:.3f}\n'.format(f_i)
    return rep_str