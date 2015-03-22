"""Code to handle the site frequency spectrum"""


class QuantizedSfs:
  """A class to store the quantized site frequency spectrum and related variant list"""
  def __init__(self, p=[float(n)/10.0 for n in range(10)], ideal_sfs=None):
    assert len(ideal_sfs) == len(p)
    self.quantized_master_list = [[] for _ in range(len(p))]
    self.quantized_p = p
    self.ideal_sfs = ideal_sfs
    assert abs(sum(self.ideal_sfs) - 1.0) < 1e3

  def add(self, variants=[]):
    """Add variants to the existing distribution, rebalancing as needed to maintain the shape of the sfs.
    :param variants: a list of lists of variants. This function does not care about the object in the list, they can be
                     None and the algorithm will still work
    """
    assert len(variants) == len(self.quantized_master_list), 'The new variants do not have the right number of buckets'
    #Do a bare add of the list
    for n in range(len(self.quantized_master_list)):
      self.quantized_master_list[n] += variants[n]
    #self.balance_sfs()  # Better to do this when we are ready

  def balance_sfs(self):
    """Pour variants from one bucket to another to make sure the sfs is approximated well"""
    variant_counts = [len(p) for p in self.quantized_master_list]
    total_variant_count = sum(variant_counts)
    ideal_variant_counts = [(f * total_variant_count) for f in self.ideal_sfs]  # sum(ideal_sfs) == 1
    # What about rounding errors?
    for n in range(len(self.quantized_master_list) - 1):
      l0 = self.quantized_master_list[n]
      if len(l0) > ideal_variant_counts[n]:
        # Easy, just pour the surplus right
        l1 = self.quantized_master_list[n + 1]
        dl = int(len(l0) - ideal_variant_counts[n] + .5)
        l1 += l0[:dl]
        self.quantized_master_list[n] = l0[dl:]
      else:
        # More complicated, need to suck surpluses from as many bins as it takes
        n2 = n + 1
        while len(l0) < ideal_variant_counts[n] and n2 < len(self.quantized_master_list):
          l1 = self.quantized_master_list[n2]
          dl = int(ideal_variant_counts[n] - len(l0) + 0.5)
          l0 += l1[:dl]
          self.quantized_master_list[n2] = l1[dl:]  # Interestingly Python handles over indexing well during slicing
          n2 += 1

  def __repr__(self):
    """Fun ASCII histogram!"""
    variant_counts = [len(p) for p in self.quantized_master_list]
    total_variant_count = sum(variant_counts)
    ideal_variant_counts = [(f * total_variant_count) for f in self.ideal_sfs]

    # We plot it as a sideways bargraph
    size_x = 80  # columns

    #Bring the data into this grid.
    scaling_factor = float(size_x) / max(variant_counts + ideal_variant_counts)
    scaled_x1 = [int(v * scaling_factor + 0.5) for v in variant_counts]
    scaled_x2 = [int(v * scaling_factor + 0.5) for v in ideal_variant_counts]
    rep_str = ''
    for v1, v2, p in zip(scaled_x1, scaled_x2, self.quantized_p):
      rep_str += '{:1.2f} '.format(p) + '-' * v1 + '\n'
      rep_str += '     ' + ' ' * v2 + '|\n'
    return rep_str


if __name__ == "__main__":
  # Couldn't resist a cool demo!
  #isfs = [n + 1 for n in range(10)]
  k = .5
  isfs = [10 * k**n for n in range(10)]
  isfs = [float(v)/sum(isfs) for v in isfs]
  sfs = QuantizedSfs(p=[float(n)/10.0 for n in range(10)], ideal_sfs=isfs)
  iv = [[n]*n*100 for n in range(10)]
  sfs.add(variants=iv)
  print(sfs)
  sfs.balance_sfs()
  print(sfs)