"""Standard population model that picks variants randomly from the master list to generate chromosomes"""
import numpy as np


class Model:
  def __init__(self, sample_size):
    """Standard population model that picks variants randomly from the master list to generate chromosomes

    :param sample_size: number of samples we should be returning
    """
    self.sample_size = sample_size
    # In more complex population models, for example simulating sexual reproduction, we would have more parameters
    # setting up things like generations to do, size of generations, number of children etc. etc.

  def samples(self, chrom_no, ml, rng_seed):
    """This returns an iterator

    :param chrom_no:  number of the chromosome being considered [1,2,3 ...]
    :param ml:        VariantList. master list of variants as created by genomes program
    :param rng_seed:  seed for random number generators
    :return: idx: A list of indexes into the variants indicating which have been chosen
    """
    rng = np.random.RandomState(rng_seed)
    gen = 0
    for n in range(self.sample_size):
      r = rng.rand(ml.variants.shape[0], 2)
      yield gen, n, ml.zip_up_chromosome(*[(r[:, 0] < ml.variants['p']).nonzero()[0], (r[:, 1] < ml.variants['p']).nonzero()[0]]), float(n) / self.sample_size
    # In more complex population models, for example simulating sexual reproduction, we would probably return an iterator
    # class that kept state representing parents etc., having worked out the population tree