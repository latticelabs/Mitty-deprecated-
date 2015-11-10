"""Code from variants.py that is a bottle neck is moved here"""
import numpy as np
cimport numpy as np
cimport cython

# (tutorial variations.json run)
#   471 ms for 10000 calls for cythonized version
# 31769 ms for 10000 calls for pure Python version
@cython.boundscheck(False)
cpdef avoid_collisions(pos, stop, idx):
  """Remove any overlapping variants from the sequence of variants indicated by idx

  :param pos:  array of start positions of master list variants
  :param stop: array of end positions of master list variants
  :param idx:  array of indexes into the variant list
  :return: an array of non-colliding indexes
  """
  cdef:
    np.ndarray[np.int32_t, ndim=1] _pos, _stop, z_idx = np.empty(len(idx), dtype=np.int32)
    np.ndarray[np.int64_t, ndim=1] _idx
    int n = 0, n2 = 0, n_max = len(idx), z_n = 0

  _pos = np.array(pos, dtype=np.int32) if type(pos) is not np.ndarray else pos
  _stop = np.array(stop, dtype=np.int32) if type(stop) is not np.ndarray else stop
  _idx = np.array(idx, dtype=np.int64) if type(idx) is not np.ndarray else idx

  while n < n_max:
    z_idx[z_n] = _idx[n]
    z_n += 1
    n2 = n + 1
    while n2 < n_max and _pos[_idx[n2]] <= _stop[_idx[n]]:
      n2 += 1  # Collision, skip
    n = n2
  return z_idx[:z_n]


# (tutorial variations.json run)
#   646 ms for 5000 calls for cythonized version
# 40090 ms for 5000 calls for pure Python version
@cython.boundscheck(False)
cpdef merge_homozygous(pos, z0, z1, bint filter_multi_allele=False):
  """Create a chromosome out of a pair of variant lists.

  :param pos:  position array from master list
  :param z0:   indexes making chrom copy 0
  :param z1:   indexes making chrom copy 1
  :param filter_multi_allele: If True discard any alleles that are both non-Ref (and not homozygous)
  :return: a list of tuples (index, genotype)
  """
  cdef:
    np.ndarray[np.int32_t, ndim=1] _pos, _z0, _z1
    int n_max0 = len(z0), n_max1 = len(z1), n0 = 0, n1 = 0, chrom_n = 0

  _pos = np.array(pos, dtype=np.int32) if type(pos) is not np.ndarray else pos
  _z0 = np.array(z0, dtype=np.int32) if type(z0) is not np.ndarray else z0
  _z1 = np.array(z1, dtype=np.int32) if type(z1) is not np.ndarray else z1

  chrom = np.empty(shape=(len(z0) + len(z1),), dtype=[('index', 'i4'), ('gt', 'i1')])
  cdef:
    np.ndarray[np.int32_t, ndim=1] chrom_index = chrom['index']
    np.ndarray[np.int8_t, ndim=1] chrom_gt = chrom['gt']

  #chrom = []
  while n0 < n_max0 and n1 < n_max1:
    if _pos[_z0[n0]] < _pos[_z1[n1]]:
      #chrom += [(_z0[n0], 0)]
      chrom_index[chrom_n], chrom_gt[chrom_n] = _z0[n0], 0
      chrom_n += 1
      n0 += 1
      continue
    if _pos[_z0[n0]] > _pos[_z1[n1]]:
      #chrom += [(_z1[n1], 1)]
      chrom_index[chrom_n], chrom_gt[chrom_n] = _z1[n1], 1
      chrom_n += 1
      n1 += 1
      continue
    # When we get here, we are equal
    if n0 < n_max0 and n1 < n_max1:  # We are equal. Are we homozygous, or just a one in a million het?
      if _z0[n0] == _z1[n1]:  # Yes, a hom
        #chrom += [(_z0[n0], 2)]
        chrom_index[chrom_n], chrom_gt[chrom_n] = _z0[n0], 2
        chrom_n += 1
      elif not filter_multi_allele:  # Just two weird hets
        #chrom += [(_z0[n0], 0)]
        #chrom += [(_z1[n1], 1)]
        chrom_index[chrom_n], chrom_gt[chrom_n] = _z0[n0], 0
        chrom_n += 1
        chrom_index[chrom_n], chrom_gt[chrom_n] = _z1[n1], 1
        chrom_n += 1
      n0 += 1
      n1 += 1  # Lets move along

  # Now zip in the remainders
  while n0 < n_max0:
    #chrom += [(_z0[n0], 0)]
    chrom_index[chrom_n], chrom_gt[chrom_n] = _z0[n0], 0
    chrom_n += 1
    n0 += 1
  while n1 < n_max1:
    #chrom += [(_z1[n1], 1)]
    chrom_index[chrom_n], chrom_gt[chrom_n] = _z1[n1], 1
    chrom_n += 1
    n1 += 1

  return chrom[:chrom_n]  #np.array((chrom_index, chrom_gt), dtype=[('index', 'i4'), ('gt', 'i1')])