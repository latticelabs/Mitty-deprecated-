import numpy as np
import h5py

import pyximport
pyximport.install(setup_args={"include_dirs": np.get_include()})
from variants_cy import *

from mitty.version import __version__


#TODO: Redesign how data is stored in hdf5 so that organization is easier to follow
# Make list of chromosomes more explicit
# Make list of samples more explicit
class Population:
  """The data is stored as follows:

  /ref_genome_meta  -> an array carrying genome metadata info from the original Fasta file
  /master_list
              /1
              /2
              ...
  /samples
          /s1
             /1
             /2
             ...
          /s2
             /1
             /2
             ...


  """
  str_dt = h5py.special_dtype(vlen=bytes)

  def __init__(self, fname='test.h5', mode='r', genome_metadata=None, in_memory=False):
    """Load a population from file, or create a new file. Over write or store the passed master list and/or samples

    :param fname:           name of the file to store/load data from.
                            If None an in-memory file is created
    :param mode:            'r'/'w' same as h5py.File modes
    :param genome_metadata: [{seq_id, seq_len, seq_md5} ...] in same order as seen in fa.gz file
                             same format as returned by Fasta.get_seq_metadata
    :param in_memory:       If True, make a file purely in memory. Mostly for testing
    """
    assert mode in ['r', 'w'], "File modes should be 'r' or 'w'"
    self.fp = h5py.File(name=fname, mode=mode,
                        driver='core' if in_memory else None, backing_store=False if in_memory else True)
    if mode not in ['r', 'r+']:  # This is an indication that we are creating a new file
      if genome_metadata is None:
        raise RuntimeError('Creating a new Population object requires genome metadata')
      self.set_genome_metadata(genome_metadata)
      self.fp.attrs['Mitty version'] = __version__

  @staticmethod
  def _ml_path(chrom):
    return '/master_list/{}'.format(chrom)

  @staticmethod
  def _s_path(sample_name=None, chrom=None):
    return '/samples/' + ((sample_name + ('/{}'.format(chrom) if chrom is not None else ''))
                          if sample_name is not None else '')

  def set_genome_metadata(self, genome_metadata):
    """Save chromosome sequence metadata

    :param genome_metadata: [{seq_id, seq_len, seq_md5} ...] in same order as seen in fa.gz file
                           same format as returned by Fasta.get_seq_metadata
    """
    dtype = [('seq_id', Population.str_dt), ('seq_len', 'i4'), ('seq_md5', Population.str_dt)]
    meta = [[gm[k] for gm in genome_metadata] for k in ['seq_id', 'seq_len', 'seq_md5']]
    self.fp.create_dataset('/ref_genome_meta', shape=(len(genome_metadata),), dtype=dtype,
                           data=np.core.records.fromarrays(meta, dtype))

  def get_genome_metadata(self):
    """Get chromosome metadata

    :returns [{seq_id, seq_len, seq_md5} ...] same format as returned by Fasta.get_seq_metadata
    """
    return [{k: x[k] for k in ['seq_id', 'seq_len', 'seq_md5']} for x in self.fp['/ref_genome_meta'][:]]

  def get_chromosome_metadata(self, chrom):
    return self.get_genome_metadata()[chrom - 1]

  def get_chromosome_list(self):
    return range(1, self.fp['/ref_genome_meta'][:].size + 1)  # Needs to be in order

  def set_master_list(self, chrom, master_list):
    """Create a new master list. Error if it already exists

    :param chrom:
    :param master_list:
    """
    assert master_list.sorted, 'Master list has not been sorted. Please check your program'
    assert len(master_list) <= 1073741823, 'Master list has more than 2^30-1 variants.'  # I want whoever gets here to mail me: kaushik.ghose@sbgenomics.com

    path = self._ml_path(chrom)
    assert path not in self.fp, "The master list exists"

    dtype = [('pos', 'i4'), ('stop', 'i4'), ('ref', Population.str_dt), ('alt', Population.str_dt), ('p', 'f2')]
    self.fp.create_dataset(name=path, shape=master_list.variants.shape,
                           dtype=dtype, data=master_list.variants, chunks=True, compression='gzip')

  def add_sample_chromosome(self, chrom, sample_name, indexes):
    """Add sample. Error if already exists

    :param chrom:  chrom number [1, 2, 3, ...]
    :param sample_name:
    :param indexes: [(chrom, gt) ...]
    """
    path = self._s_path(sample_name, chrom)
    assert path not in self.fp, "This sample/chrom exists"
    assert self._ml_path(chrom) in self.fp, "This chromosome is absent in the master list"

    self.fp.create_dataset(name=path, shape=indexes.shape, dtype=[('index', 'i4'), ('gt', 'i1')],
                           data=indexes, chunks=True, compression='gzip')

  def get_variant_master_list_count(self, chrom):
    path = self._ml_path(chrom)
    return self.fp[path].size if path in self.fp else 0

  def get_variant_master_list(self, chrom):
    """Return the whole master variant list for this chromosome"""
    ml = VariantList()
    if self._ml_path(chrom) in self.fp:
      ml.variants = self.fp[self._ml_path(chrom)][:]
    return ml

  def get_sample_variant_count(self, chrom, sample_name):
    path = self._s_path(sample_name, chrom)
    return self.fp[path].size if path in self.fp else 0

  def get_sample_variant_index_for_chromosome(self, chrom, sample_name):
    """Return the indexes pointing to the master list for given sample and chromosome"""
    path = self._s_path(sample_name, chrom)
    return self.fp[path][:] if path in self.fp else np.array([], dtype=[('index', 'i4'), ('gt', 'i1')])

  def get_sample_variant_list_for_chromosome(self, chrom, sample_name, ignore_zygosity=False):
    """Return variant list for this sample and chromosome."""
    ml = self.get_variant_master_list(chrom)
    v_idx = self.get_sample_variant_index_for_chromosome(chrom, sample_name)
    if ignore_zygosity:
      return ml.variants[v_idx['index']]
    else:
      return [ml.variants[v_idx['index'][(v_idx['gt'] == 0) | (v_idx['gt'] == 2)]],
              ml.variants[v_idx['index'][(v_idx['gt'] == 1) | (v_idx['gt'] == 2)]]]

  def get_sample_names(self):
    """Return a list of sample names"""
    return self.fp[self._s_path()].keys()

  def get_version(self):
    return self.fp.attrs['Mitty version']

  #TODO: make more detailed
  def __repr__(self):
    """Pretty print the genome file"""
    return self.pretty_print_summary()

  def pretty_print_summary(self, sample_name=None):
    """Give us a nice printed representation of the population file"""
    rep_str = """
    ---------------------------------------
    Genome file. Mitty version {mv:s}
    ---------------------------------------
    {chrom_cnt:d} chromosomes
    {sample_cnt:d} samples\n
    """.format(mv=self.get_version(), chrom_cnt=len(self.get_chromosome_list()), sample_cnt=len(self.get_sample_names()))
    pop_v_cnt = [self.get_variant_master_list_count(chrom=chrom) for chrom in self.get_chromosome_list()]
    if sample_name is not None:
      sample_v_cnt = [self.get_sample_variant_count(chrom=chrom, sample_name=sample_name) for chrom in self.get_chromosome_list()]
    else:
      sample_v_cnt = []

    # Now print it out nicely
    rep_str += 'Variant counts' + (' (sample: {:s})\n'.format(sample_name) if sample_name is not None else '\n')
    rep_str += '\tChrom\tPop' + ('\t\tSample\n' if sample_name is not None else '\n')
    rep_str += '\t     \tCount' + ('\t\tCount\n' if sample_name is not None else '\n')
    rep_str += '\t-----\t-----' + ('\t\t-----\n' if sample_name is not None else '\n')
    for n, chrom in enumerate(self.get_chromosome_list()):
      rep_str += '\t{:d}\t{:d}'.format(chrom, pop_v_cnt[n]) + ('\t\t{:d}\n'.format(sample_v_cnt[n]) if sample_name is not None else '\n')
    rep_str += '\t-----\t-----' + ('\t\t-----\n' if sample_name is not None else '\n')
    rep_str += '\tTotal\t{:d}'.format(sum(pop_v_cnt)) + ('\t\t{:d}\n'.format(sum(sample_v_cnt)) if sample_name is not None else '\n')
    return rep_str


def l2ca(l):
  """Convenience function that converts a Python list of tuples into an numpy structured array corresponding to a
  chromosome index array"""
  return np.array(l, dtype=[('index', 'i4'), ('gt', 'i1')])


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

  def zip_up_chromosome(self, idx0, idx1, filter_multi_allele=False):
    """Given two chromosomes, go through each copy, variant by variant, making sure they don't clash and merging any
    homozygous ones. Return us a chromosome array. Will sort master list if not sorted
    :param idx0: proposed indexes for chrom copy 0
    :param idx1: proposed indexes for chrom copy 1
    :param filter_multi_allele: If True discard any alleles that are both non-Ref (and not homozygous)
    :returns: chrom, a list of tuples (index, genotype)
    """
    if not self.sorted:
      self.sort()

    # Pass 1: get rid of the colliding variants
    pos, stop = self.variants['pos'], self.variants['stop']
    z0, z1 = avoid_collisions(pos, stop, idx0), avoid_collisions(pos, stop, idx1)

    # Pass 2: merge homozygous where needed
    return merge_homozygous(pos, z0, z1, filter_multi_allele)

  def generate_chromosome(self, rng):
    """Convenient wrapper around select and zip_up_chromosome
    :param rng: a random number generator with the .rand method returning uniform random numbers (0.0, 1.0)
    :returns: chrom, a list of tuples (index, genotype)
    """
    return self.zip_up_chromosome(*self.select(rng))

  def __repr__(self):
    """Fun ASCII histogram!"""
    if self.variants.shape[0] == 0:
      return '<empty>'

    if self.site_freq_spectrum is not None:
      sfs_p, sfs = self.site_freq_spectrum
      ideal_cnt = [f * self.variants.shape[0] for f in sfs]
    else:
      sfs_p = np.linspace(0, self.variants['p'].max(), num=11)  # Default is to histogram in 11 bins
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


# This is the pure python version. The cythonized version is in variants_cy.pyx
def py_avoid_collisions(pos, stop, idx):
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


# This is the pure python version. The cythonized version is in variants_cy.pyx
def py_merge_homozygous(pos, z0, z1, filter_multi_allele=False):
  """Create a chromosome out of a pair of variant lists.

  :param pos:  position array from master list
  :param z0:   indexes making chrom copy 0
  :param z1:   indexes making chrom copy 1
  :param filter_multi_allele: If True discard any alleles that are both non-Ref (and not homozygous)
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
      elif not filter_multi_allele:  # Just two weird hets
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

  return np.array(chrom, dtype=[('index', 'i4'), ('gt', 'i1')])