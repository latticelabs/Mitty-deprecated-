"""Standard population model that picks variants randomly from the master list based on probability value"""
import numpy as np
import json

__example_param_text = """
{
  "standard": {
    "sample_size": 10,
    "force_homozygous": True, # Force sample variants to be 1|1
    "filter_multi_allele": False,  # Take out locii with different variants on the two copies
    "filter_hom": False,  # Take out homozygous
    "max_v_count": 100,  # Maximum number of variants
    "min_v_spacing": 1000  # Minimum gap between variants on same copy
  }
}
"""

_description = __doc__ + '\nExample parameters:\n' + __example_param_text

#_example_params = json.loads(__example_param_text)
_example_params = eval(__example_param_text)


class Model:
  def __init__(self, sample_size=10, force_homozygous=False, filter_multi_allele=False, filter_hom=False,
               max_v_count=None, min_v_spacing=None):
    """Standard population model that picks variants randomly from the master list to generate chromosomes

    :param sample_size: number of samples we should be returning
    :param force_homozygous: force all variants to be homozygous
    """
    self.sample_size = sample_size
    self.force_homozygous = force_homozygous
    self.filter_multi_allele = filter_multi_allele
    self.filter_hom = filter_hom
    self.max_v_count = max_v_count
    self.min_v_spacing = min_v_spacing
    # In more complex population models, for example simulating sexual reproduction, we would have more parameters
    # setting up things like generations to do, size of generations, number of children etc. etc.

  def filter_sample(self, r, ml=None, rng=None):
    """Apply ad hoc filtering to the generated samples to fit criteria tester wants"""
    # The following filtering operations may be slow!
    if self.force_homozygous:
      r[:, 1] = r[:, 0]
    if self.max_v_count is not None:
      all_idx = ((r[:, 0] < ml.variants['p']) | (r[:, 1] < ml.variants['p'])).nonzero()[0]
      if all_idx.size > self.max_v_count:
        r[rng.choice(all_idx, size=(all_idx.size - self.max_v_count), replace=False), :] = 1.0
    if self.filter_hom:
      raise NotImplementedError
    if self.min_v_spacing:  # Expensive
      starts = ml.variants['pos']
      stops = ml.variants['stop']
      for cpy in [0, 1]:
        idx = (r[:, cpy] < ml.variants['p']).nonzero()[0]
        n0, n1 = 0, 1
        while n0 < n1 < idx.size - 1:
          while starts[idx[n1]] - stops[idx[n0]] < self.min_v_spacing:
            r[idx[n1], cpy] = 1.0  # Puts this out of contention for selection
            if n1 < idx.size - 1:
              n1 += 1
            else:
              break
          n0 = n1
          n1 += 1
    return r

  def samples(self, chrom_no=None, ml=None, rng_seed=1):
    """This returns an iterator

    :param chrom_no:  number of the chromosome being considered [1,2,3 ...]  (ignored here)
    :param ml:        VariantList. master list of variants as created by genomes program
    :param rng_seed:  seed for random number generators
    :return: A generator returning (generation no, serial_no, chromosome, % samples done) for each sample in population
    """
    rng = np.random.RandomState(rng_seed)
    if self.force_homozygous or self.filter_multi_allele or self.filter_hom or \
      self.max_v_count is not None or self.min_v_spacing is not None:
      ad_hoc_filtering = True
    else:
      ad_hoc_filtering = False

    gen = 0
    for n in range(self.sample_size):
      r = rng.rand(ml.variants.shape[0], 2)
      if ad_hoc_filtering:
        r = self.filter_sample(r, ml, rng)
      yield 'g{:d}_s{:d}'.format(gen, n), ml.zip_up_chromosome(*[(r[:, 0] < ml.variants['p']).nonzero()[0], (r[:, 1] < ml.variants['p']).nonzero()[0]], filter_multi_allele=self.filter_multi_allele), float(n + 1) / self.sample_size
    # In more complex population models, for example simulating sexual reproduction, we would probably return an iterator
    # class that kept state representing parents etc., having worked out the population tree

  def get_sample_count_estimate(self):
    """Give us an as exact as possible estimate of how many samples we will produce"""
    return self.sample_size