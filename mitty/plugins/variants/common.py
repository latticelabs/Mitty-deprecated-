import numpy as np


def scale_probability_and_validate(self_p, p, f):
  """Scale the per base per sample probability to a per base per population value if the site frequency spectrum is
  given. Raise all sorts of assertion errors if things go wonky.

  :param self_p:   the per base, per sample probability value
  :param p:        the site frequency spectrum values
  :param f:
  :return: p_eff  - the effective (per population) probability
  """
  if p is not None:
    assert type(p) == np.ndarray, 'Site Freq spectrum should be numpy array'
    assert type(f) == np.ndarray, 'Site Freq spectrum should be numpy array'
    assert abs(1.0 - f.sum()) < 1e-6, 'Site freq spectrum should sum to 1.0'
    p_eff = self_p / (p * f).sum()  # Effective per-base probability for master list
    assert 0 <= p_eff <= 1.0, 'We are getting an illegal effective probability value. Recheck site freq spectrum and p and see manual'
  else:
    p_eff = self_p
  return p_eff