"""Standard population model that picks variants randomly from the master list based on probability value"""
import numpy as np
import json

__example_param_text = """
{
  "standard": {
    "sample_size": 10,
    "force_homozygous": true
  }
}
"""

_description = __doc__ + '\nExample parameters:\n' + __example_param_text

_example_params = json.loads(__example_param_text)


class Model:
  def __init__(self, sample_size=10, force_homozygous=False):
    """Standard population model that picks variants randomly from the master list to generate chromosomes

    :param sample_size: number of samples we should be returning
    :param force_homozygous: force all variants to be homozygous
    """
    self.sample_size = sample_size
    self.force_homozygous = force_homozygous
    # In more complex population models, for example simulating sexual reproduction, we would have more parameters
    # setting up things like generations to do, size of generations, number of children etc. etc.

  def samples(self, chrom_no=None, ml=None, rng_seed=1):
    """This returns an iterator

    :param chrom_no:  number of the chromosome being considered [1,2,3 ...]  (ignored here)
    :param ml:        VariantList. master list of variants as created by genomes program
    :param rng_seed:  seed for random number generators
    :return: A generator returning (generation no, serial_no, chromosome, % samples done) for each sample in population
    """
    rng = np.random.RandomState(rng_seed)
    gen = 0
    for n in range(self.sample_size):
      if self.force_homozygous:
        r = np.tile(rng.rand(ml.variants.shape[0], 1), 2)
      else:
        r = rng.rand(ml.variants.shape[0], 2)
      yield 'g{:d}_s{:d}'.format(gen, n), ml.zip_up_chromosome(*[(r[:, 0] < ml.variants['p']).nonzero()[0], (r[:, 1] < ml.variants['p']).nonzero()[0]]), float(n + 1) / self.sample_size
    # In more complex population models, for example simulating sexual reproduction, we would probably return an iterator
    # class that kept state representing parents etc., having worked out the population tree

  def get_sample_count_estimate(self):
    """Give us an as exact as possible estimate of how many samples we will produce"""
    return self.sample_size