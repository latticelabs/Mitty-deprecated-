"""A population model that creates samples with more and more variants. Suitable for the aligner paper experiments
^ = intersection
E = subset

vx ^ v0 = v0
vx ^ v1 = v0
...
vx ^ vn = v0

v0 E v1
v1 E v2
v2 E v3
...
v(n-1) E vn

This plugin does not honor the site frequency spectrum model and ignores the original 'p' values
"""
import numpy as np

__example_param_text = """
{
  "standard": {
    "p_vx": 0.2,
    "p_vn": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
  }
}
"""

_description = __doc__ + '\nExample parameters:\n' + __example_param_text

#_example_params = json.loads(__example_param_text)
_example_params = eval(__example_param_text)


class Model:
  def __init__(self, p_vx, p_vn):
    """A population model that creates samples with more and more variants. Suitable for the aligner paper experiments

    :param p_vx: probability value for vx set
    :param p_vn: probability values for v0, v1, v2, v3 .... set
    """
    self.p_vx, self.p_vn = p_vx, p_vn

  def samples(self, chrom_no=None, ml=None, rng_seed=1, **kwargs):
    """This returns an iterator

    :param chrom_no:  number of the chromosome being considered [1,2,3 ...]  (ignored here)
    :param ml:        VariantList. master list of variants as created by genomes program
    :param rng_seed:  seed for random number generators
    :return: A generator returning (generation no, serial_no, chromosome, % samples done) for each sample in population

    Algorithm: (Repeat for each chromosome copy)

    Generate random numbers r same size as variants list
    Select vx <= r < p_vx
    Pick a random subset of v0 as v1 size(v1)/size(v0) = p_v1/p_v0
    Set all r corresponding to v0 - v1 as 1.0 so we never select these again
    Pick v2, v3 ... by comparing r to p_v2, p_v3 and so on


    """
    assert 0 <= self.p_vx <= 1.0, 'p_vx needs to be >= 0 and <= 1.0'
    assert self.p_vx > self.p_vn[0], 'p_vx needs to be > p_vn[0]'
    for n in range(len(self.p_vn) - 1):
      assert self.p_vn[n] < self.p_vn[n + 1], 'p_vn needs to be in ascending order'
      assert 0 <= self.p_vn[n] <= 1.0, 'p_vn needs to be >= 0 and <= 1.0'

    rng = np.random.RandomState(rng_seed)
    r = rng.rand(ml.variants.shape[0], 2)

    idx_vx = [None, None]
    for cpy in [0, 1]:
      idx_vx[cpy] = np.sort(rng.choice(ml.variants.shape[0], size=int(ml.variants.shape[0] * self.p_vx), replace=False))
      # Take elements in vx that are not going to be in v0 completely out of circulation
      r[idx_vx[cpy][(r[idx_vx[cpy]] >= self.p_vn[0]).nonzero()[0]], cpy] = 1.1
      # Now all elements for r < 1.0 are either in vx ^ v0 or not in vx

    for n in range(len(self.p_vn) + 1):
      if n == 0:
        this_idx, sample_name = idx_vx, 'vx'
      else:
        this_idx, sample_name = [(r[:, cpy] < self.p_vn[n - 1]).nonzero()[0] for cpy in [0, 1]], 'v{:d}'.format(n - 1)

      yield sample_name, ml.zip_up_chromosome(*this_idx), float(n + 1) / self.get_sample_count_estimate()

  def get_sample_count_estimate(self):
    """Give us an as exact as possible estimate of how many samples we will produce"""
    return 1 + len(self.p_vn)