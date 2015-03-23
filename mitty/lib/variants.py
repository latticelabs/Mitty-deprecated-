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

  def balance_probabilities(self, p, sfs):
    """Use the ideal site probability spectrum to rescale the probability values
    :param p: probability values
    :param sfs: proportion"""
    raise NotImplementedError('To do')
    #idx = self.variants['p'].argsort()  # We need the data sorted by probability value for this to work

  def select(self, rng):
    """
    :param rng: a random number generator with the .rand method returning uniform random numbers (0.0, 1.0)
    :return: idx: A list of indexes into the variants indicating which have been chosen
    """
    r = rng.rand(self.variants.shape[0], 2)
    return [(r[:, 0] < self.variants['p']).nonzero()[0], (r[:, 1] < self.variants['p']).nonzero()[0]]

  def zip_up_chromosome(self, idx0, idx1):
    """Given two chromosomes, go through each copy, variant by variant, making sure they don't clash and merging any
    homozygous ones. Return us a chromosome array. This should only be done if the list is sorted"""
    if not self.sorted:
      self.sort()

    # Pass 1: get rid of the colliding variants
    pos, stop = self.variants['pos'], self.variants['stop']
    z0, z1 = avoid_collisions(pos, stop, idx0), avoid_collisions(pos, stop, idx1)

    # Pass 2: merge homozygous where needed
    return merge_homozygous(pos, z0, z1)

  def generate_chromosome(self, rng):
    """Wrapper around select and zip_up_chromosome."""
    return self.zip_up_chromosome(*self.select(rng))


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
    while n0 < n_max0 and n1 < n_max1 and pos[z0[n0]] < pos[z1[n1]]:
      chrom += [(z0[n0], 0)]
      n0 += 1
    while n0 < n_max0 and n1 < n_max1 and pos[z1[n1]] < pos[z0[n0]]:
      chrom += [(z1[n1], 1)]
      n1 += 1
    # When we get here, we are either equal, or we've run out of stuff
    if n0 < n_max0 and n1 < n_max1:  # We are equal. Are we homozygous, or just a one in a million het?
      if z0[n0] == z1[n1]:  # Yes, a hom
        chrom += [(z0[n0], 2)]
      else:  # Just two weird hets
        chrom += [(z0[n0], 0)]
        chrom += [(z1[n1], 1)]
      n0 += 1
      n1 += 1  # Lets move along
  return chrom





