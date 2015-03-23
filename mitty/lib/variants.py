import numpy as np


class VariantList:
  """Use a numpy recarray to store a list of variants."""
  def __init__(self, pos_a=[], stop_a=[], ref_a=[], alt_a=[], p_a=[]):
    #self.variant = np.zeros((len(pos_a),), dtype=[('pos', 'i4'), ('stop', 'i4'), ('ref', 'object'), ('alt', 'object'), ('p', 'f2')])
    #self.variant['pos'], self.variant['stop'], self.variant['ref'], self.variant['alt'], self.variant['p'] = pos_a, stop_a, ref_a, alt_a, p_a
    #self.pos, self.stop, self.ref, self.alt, self.p = pos_a, stop_a, ref_a, alt_a, p_a
    self.variants = np.core.records.fromarrays([pos_a, stop_a, ref_a, alt_a, p_a], dtype=[('pos', 'i4'), ('stop', 'i4'), ('ref', 'object'), ('alt', 'object'), ('p', 'f2')])

  def add(self, pos_a=[], stop_a=[], ref_a=[], alt_a=[], p_a=[]):
    new_variants = np.core.records.fromarrays([pos_a, stop_a, ref_a, alt_a, p_a], dtype=[('pos', 'i4'), ('stop', 'i4'), ('ref', 'object'), ('alt', 'object'), ('p', 'f2')])
    self.variants = np.concatenate((self.variants, new_variants))

  def sort(self):
    """Sort us in place by the position."""
    idx = self.variants['pos'].argsort()  # recarray sort uses all fields and is wasteful
    self.variants = self.variants[idx]
    #self.pos, self.stop, self.ref, self.alt, self.p = self.pos[idx], self.stop[idx], self.ref[idx], self.alt[idx], self.p[idx]

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
    homozygous ones. Return us a chromosome array"""
    def avoid_collisions(pos, stop, idx):
      n, n_max = 0, len(idx)
      z_idx = []
      while n < n_max - 1:
        z_idx += [idx[n]]
        if pos[idx[n + 1]] < stop[idx[n]]:
          n += 2  # Collision, skip
        else:
          n += 1  # Add next
      return z_idx

    # Pass 1: get rid of the colliding variants
    pos, stop = self.variants['pos'], self.variants['stop']
    z0, z1 = avoid_collisions(pos, stop, idx0), avoid_collisions(pos, stop, idx1)

    # Merge homozygous
    n_max0, n_max1 = len(z0), len(z1)
    n0, n1 = 0, 0
    chrom = []
    while n0 < n_max0 and n1 < n_max1:
      while pos[z0[n0]] < pos[z1[n1]] and n0 < n_max0 and n1 < n_max1:
        chrom += [(z0[n0], 0)]
        n0 += 1
      while pos[z1[n1]] < pos[z0[n0]] and n0 < n_max0 and n1 < n_max1:
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

  def generate_chromosome(self, rng):
    """Wrapper around select and zip_up_chromosome."""
    return self.zip_up_chromosome(*self.select(rng))







