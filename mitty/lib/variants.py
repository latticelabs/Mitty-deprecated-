import numpy as np


class VariantList:
  """Use a numpy recarray to store a list of variants."""
  def __init__(self, pos_a=[], stop_a=[], ref_a=[], alt_a=[], p_a=[]):
    """Initialize the array from individual arrays/lists

    :param pos_a: position vector
    :param stop_a: stop vector
    :param ref_a: reference bases
    :param alt_a: alt bases
    :param p_a: probability value for the variants
    :return: Object
    """
    self.variants = np.core.records.fromarrays([pos_a, stop_a, ref_a, alt_a, p_a],
                                               dtype=[('pos', 'i4'), ('stop', 'i4'), ('ref', 'object'), ('alt', 'object'), ('p', 'f2')])
    self.sorted = False
    self.site_freq_spectrum = None

  def __len__(self):
    return self.variants.shape[0]

  def add(self, pos_a=[], stop_a=[], ref_a=[], alt_a=[], p_a=[]):
    """Add more variants to the list

    :param pos_a: position vector
    :param stop_a: stop vector
    :param ref_a: reference bases
    :param alt_a: alt bases
    :param p_a: probability value for the variants
    :return: Object
    """
    new_variants = np.core.records.fromarrays([pos_a, stop_a, ref_a, alt_a, p_a],
                                              dtype=[('pos', 'i4'), ('stop', 'i4'), ('ref', 'object'), ('alt', 'object'), ('p', 'f2')])
    self.variants = np.concatenate((self.variants, new_variants))
    self.sorted = False

  def sort(self):
    """Sort us in place by the position."""
    idx = self.variants['pos'].argsort()  # recarray sort uses all fields and is wasteful
    self.variants = self.variants[idx]
    self.sorted = True

  def balance_probabilities(self, p, f):
    """Use the ideal site probability spectrum to rescale the probability values
    :param p: probability values
    :param f: proportion

    sum(f) = 1.0 for this to work"""
    assert len(p) == len(f)
    assert abs(1.0 - sum(f)) < 1e-6
    idx = self.variants['p'].argsort()  # We need the data sorted by probability value for this to work
    n_max = len(self)
    n = 0
    for p_i, f_i in zip(p, f):
      self.variants['p'][idx[n:n + int(f_i * n_max + .5)]] = p_i  # Over index is handled gracefully
      n += int(f_i * n_max + .5)
    self.site_freq_spectrum = (p, f)

  def select(self, rng):
    """Use the rng to select variants from the master list based on their probabilities
    :param rng: a random number generator with the .rand method returning uniform random numbers (0.0, 1.0)
    :return: idx: A list of indexes into the variants indicating which have been chosen
    """
    r = rng.rand(self.variants.shape[0], 2)
    return [(r[:, 0] < self.variants['p']).nonzero()[0], (r[:, 1] < self.variants['p']).nonzero()[0]]

  def zip_up_chromosome(self, idx0, idx1):
    """Given two chromosomes, go through each copy, variant by variant, making sure they don't clash and merging any
    homozygous ones. Return us a chromosome array. Will sort master list if not sorted
    :param idx0: proposed indexes for chrom copy 0
    :param idx1: proposed indexes for chrom copy 1
    :returns: chrom, a list of tuples (index, genotype)
    """
    if not self.sorted:
      self.sort()

    # Pass 1: get rid of the colliding variants
    pos, stop = self.variants['pos'], self.variants['stop']
    z0, z1 = avoid_collisions(pos, stop, idx0), avoid_collisions(pos, stop, idx1)

    # Pass 2: merge homozygous where needed
    return merge_homozygous(pos, z0, z1)

  def generate_chromosome(self, rng):
    """Convenient wrapper around select and zip_up_chromosome
    :param rng: a random number generator with the .rand method returning uniform random numbers (0.0, 1.0)
    :returns: chrom, a list of tuples (index, genotype)
    """
    return self.zip_up_chromosome(*self.select(rng))

  def __repr__(self):
    """Fun ASCII histogram!"""
    if self.site_freq_spectrum is not None:
      sfs_p, sfs = self.site_freq_spectrum
      ideal_cnt = [f * self.variants.shape[0] for f in sfs]
    else:
      sfs_p = np.linspace(0, 1.0, num=11)  # Default is to histogram in 11 bins
      ideal_cnt = [0 for _ in range(11)]

    # Now histogram the actual data
    dp = (sfs_p[1:] - sfs_p[:-1]) / 2.0 if len(sfs_p) > 1 else [0.5]
    actual_cnt, be = np.histogram(self.variants['p'], np.concatenate(([0], sfs_p[:-1] + dp, [sfs_p[-1] + dp[-1]])))

    # We plot it as a sideways bar-graph
    size_x = min(max(actual_cnt), 80)  # columns

    #Bring the data into this grid.
    scaling_factor = float(size_x) / max(max(actual_cnt), max(ideal_cnt))
    scaled_actual = [int(v * scaling_factor + 0.5) for v in actual_cnt]
    scaled_ideal = [int(v * scaling_factor + 0.5) for v in ideal_cnt]
    rep_str = ''
    for na, sc_a, sc_i, p in zip(actual_cnt, scaled_actual, scaled_ideal, sfs_p):
      rep_str += '{:1.2f} '.format(p)
      if sc_i <= sc_a:  # The | for the ideal comes before or overlaps with the last -
        rep_str += '-' * (sc_i - 1) + ('|' if sc_i else '') + '-' * (sc_a - sc_i) + ' {:d}\n'.format(na)
      else:  # The | comes beyond the last -
        rep_str += '-' * sc_a + ' ' * (sc_i - sc_a) + '| {:d}\n'.format(na)
    return rep_str


def avoid_collisions(pos, stop, idx):
  """Remove any overlapping variants from the sequence of variants indicated by idx

  :param pos:  array of start positions of master list variants
  :param stop: array of end positions of master list variants
  :param idx:  array of indexes into the variant list
  :return: an array of non-colliding indexes
  """
  n, n_max = 0, len(idx)
  z_idx = []
  while n < n_max:
    z_idx += [idx[n]]
    n2 = n + 1
    while n2 < n_max and pos[idx[n2]] <= stop[idx[n]]:
      n2 += 1  # Collision, skip
    n = n2
  return z_idx


def merge_homozygous(pos, z0, z1):
  """Create a chromosome out of a pair of variant lists.

  :param pos:  position array from master list
  :param z0:   indexes making chrom copy 0
  :param z1:   indexes making chrom copy 1
  :return: a list of tuples (index, genotype)
  """
  n_max0, n_max1 = len(z0), len(z1)
  n0, n1 = 0, 0
  chrom = []
  while n0 < n_max0 and n1 < n_max1:
    if pos[z0[n0]] < pos[z1[n1]]:
      chrom += [(z0[n0], 0)]
      n0 += 1
      continue
    if pos[z0[n0]] > pos[z1[n1]]:
      chrom += [(z1[n1], 1)]
      n1 += 1
      continue
    # When we get here, we are equal
    if n0 < n_max0 and n1 < n_max1:  # We are equal. Are we homozygous, or just a one in a million het?
      if z0[n0] == z1[n1]:  # Yes, a hom
        chrom += [(z0[n0], 2)]
      else:  # Just two weird hets
        chrom += [(z0[n0], 0)]
        chrom += [(z1[n1], 1)]
      n0 += 1
      n1 += 1  # Lets move along

  # Now zip in the remainders
  while n0 < n_max0:
    chrom += [(z0[n0], 0)]
    n0 += 1
  while n1 < n_max1:
    chrom += [(z1[n1], 1)]
    n1 += 1

  return chrom