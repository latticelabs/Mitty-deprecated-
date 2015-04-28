import pyximport; pyximport.install()
from util_cython import *


def growth_curve_from_sfs(p, f, n_list=[1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000]):
  """Given the site frequency spectrum, compute the growth curve of the variant count

  :param p: array of probability values
  :param f: array of site frequencies
  :param max_samples: number of samples to compute the curve to
  :return: growth curve with max_samples elements going from 0 to 1
  """
  return n_list, [(f * (1.0 - (1.0 - p) ** n)).sum() for n in n_list]