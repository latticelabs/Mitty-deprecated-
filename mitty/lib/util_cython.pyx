import string
from itertools import izip

import cython
import numpy as np
cimport numpy as np
cimport cython

from mitty.lib import SEED_MAX


DNA_complement = string.maketrans('ATCGN', 'TAGCN')


def initialize_rngs(unsigned long master_seed, int n_rngs=4):
  """Return n_rngs initialized from the master_seed"""
  return [np.random.RandomState(seed=seed)
          for seed in np.random.RandomState(seed=master_seed).randint(SEED_MAX, size=n_rngs)]


def place_poisson_seq(rng, float p, unsigned long start_x, unsigned long end_x, bytes seq):
  """Given a random number generator, a probability and an end point, generate poisson distributed events. Skip bases
  that are 'N'.  For short end_p this may, by chance, generate fewer locations that normal"""
  if p == 0.0:
    return np.array([])

  cdef:
    char *s = seq
    unsigned long est_block_size = <unsigned long>(<float>end_x * p * 1.2)
    unsigned long idx

  these_locs = rng.geometric(p=p, size=est_block_size).cumsum()
  return np.array([idx for idx in these_locs[np.searchsorted(these_locs, start_x):np.searchsorted(these_locs, end_x)] if s[idx] != 'N'], dtype='i4')


def discard_deletions_in_illegal_regions(
    bytes ref, np.ndarray[np.int32_t, ndim=1] start_loc, np.ndarray[np.int32_t, ndim=1] stop_loc):
  """Return lists of del_locs, del_ends, refs, alts for deletes that ensure
     1. The deletions don't go past the end of ref
     2. There are no 'N's anywhere. The FASTA reader converts all IUPAC codes to 'N's"""
  cdef:
    char *s = ref
    np.ndarray[np.int32_t, ndim=1] idx = np.empty(start_loc.size, dtype=np.int32)
    int n = 0, n_idx = 0, n_max = start_loc.size, ref_len = len(ref)

  while n < n_max:
    if legal_del(s, start_loc[n], stop_loc[n], ref_len):
      idx[n_idx] = n
      n_idx += 1
    n += 1

  if n_idx > 0:
    _start, _stop = start_loc[idx[:n_idx]], stop_loc[idx[:n_idx]]
    refs, alts = map(list, izip(*((ref[st:nd], ref[st]) for (st, nd) in np.nditer([_start, _stop]))))
  else:
    _start, _stop, refs, alts = np.array([], dtype='i4'), np.array([], dtype='i4'), [], []

  return _start, _stop, refs, alts


cdef int legal_del(char *s, int start, int stop, int s_len):
  """Must be within length of s, can't contain an 'N'"""
  cdef:
    int n = start
  if stop > s_len:
    return 0

  while n < stop:
    if s[n] == 'N':
      return 0
    n += 1

  return 1




cdef unsigned char sub_base(unsigned char orig_base, unsigned char sub_mat[85][3], float ct_mat[85][3], float r):
  # sub_matrix
  #          0  1  2
  # A (65)   C  G  T
  # C (67)   A  G  T
  # G (71)   A  C  T
  # T (84)   A  C  G

  cdef:
    unsigned char i
  for i in range(3):
    if r < ct_mat[orig_base][i]:
      break
  return sub_mat[orig_base][i]


def base_subs(bytes seq, sub_pts, t_mat, rng):
  """
  t_mat ->
     A C G T
  A  . x x x
  C  x . x x
  G  x x . x
  T  x x x .
  """
  cdef:
    char *s = seq
    unsigned char i, j, ob, sb
    unsigned char sub_mat[85][3]
    float ct_mat[85][3]

  # Create the ct_mat from the basic t_mat
  for i, ob in enumerate([65, 67, 71, 84]):
    jj = 0
    for j, sb in enumerate([65, 67, 71, 84]):
      if i == j: continue
      ct_mat[ob][jj] = t_mat[i][j]  # First pass, compact the matrix
      sub_mat[ob][jj] = sb  #
      jj += 1
    # Second pass, do cumulative probabilities
    ct_mat[ob][1] += ct_mat[ob][0]
    ct_mat[ob][2] += ct_mat[ob][1]

  r = rng.rand(len(sub_pts))
  cdef unsigned long q
  return [<bytes> sub_base(s[sub_pts[q]], sub_mat, ct_mat, r[q]) for q in range(len(sub_pts))]


def add_p_end_to_t_mat(t_mat, p_end):
  """Given p_end, incorporate it into the t_mat."""
  p_end_1 = 1.0 - p_end
  return [[p_end_1 * (t_mat[i][j]/sum(t_mat[i])) for j in range(4)] + [p_end] for i in range(4)]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef markov_chain_sequence_gen(
    char first_letter, double ct_mat[4][5], unsigned char *seq, unsigned long int *l, unsigned long int max_len, rng):
  """sequence_gen(float t_mat[4][5], char *alphabet[4])
  :param (char) first_letter: the first letter of the sequence
  :param (float) ct_mat: 4x5 cumulative transition probability matrix (ACGTx)
  :param (char*) seq: a max_len long string allocation. Filled out with sequence
  :param (int*) l: length of sequence
  :param (int) max_len: We pinch off a sequence that is this long
  :param rng: numpy random number generator object that has rand

  Cumulative Transition Probability matrix (x represents the break state - end of string)

     A C G T x
  A  . . . . .
  C  . . . . .
  G  . . . . .
  T  . . . . .

  Each row indicates the cumulative probability of character P (row) going to character Q (column) and thus ends with 1.
  We use cumulative probability rather than actual probability to make our code faster - we just check if the random
  number we generated is less than col0, col1, col2 and so on...
  """
  cdef:
    unsigned char last_letter = 0, n
    const char *alphabet = "ACGT"
    np.ndarray[double, ndim=1, mode="c"] rnd = rng.rand(max_len + 5)  # Give a little overage for premature exit states
    Py_ssize_t rnd_idx = 0
    double r
    bint keep_running = 1

  for n in range(4):
    if first_letter == alphabet[n]:
      last_letter = n
      break

  seq[0] = first_letter
  l[0] = 1

  while keep_running:
    if rnd_idx == rnd.shape[0]:  # We're out of randomness. Need more. Need more.
      rnd = rng.rand(max_len)
      rnd_idx = 0
    r = rnd[rnd_idx]
    rnd_idx += 1
    for n in range(5):
      if r <= ct_mat[last_letter][n]:
        break
    if n == 4:
      if l[0] > 1: keep_running = 0 # We can stop if we have a length 2 sequence
    else:
      seq[l[0]] = alphabet[n]
      last_letter = n
      l[0] += 1

    if l[0] == max_len + 1: keep_running = 0


def markov_sequences(bytes seq, ins_pts, max_len, t_mat, rng):
  """markov_sequences(seq, ins_pts, max_len, t_mat, rng)
  Return us insertions at the requested positions.
  :param (str) seq: the reference sequence. Needed for first letters of insertions
  :param ins_pts: list/array of insertion points
  :param max_len: either a scalar or list/array (same length as ins_pts) that
  :param (4x5 list) t_mat: transition matrix, including prob of termination
  :param rng: numpy random number generator object that has rand
  :returns a list of insertion sequences
  """
  if type(max_len) is int:
    max_max_len = max_len
    max_lens = [max_len] * len(ins_pts)
  else:
    max_lens = max_len
    assert len(max_len) == len(ins_pts), 'Lengths of insertion points and max lengths must be equal'
    max_max_len = max(max_len) if len(max_len) else 1

  pre_alloc_str = 'N' * (max_max_len + 1)  # Pre-allocated string. +1 to include ref base

  cdef:
    unsigned char *s = seq
    double ct_mat[4][5]
    unsigned long int l
    unsigned char *pre_string = pre_alloc_str

  # Convert transition matrix to cumulative transition matrix
  for i in range(4):
    ct_mat[i][0] = t_mat[i][0]
    for j in range(1, 5):
      ct_mat[i][j] = ct_mat[i][j-1] + t_mat[i][j]

  insertions, lengths = [], []
  for ip, max_len in izip(ins_pts, max_lens):
    markov_chain_sequence_gen(s[ip], ct_mat, pre_string, &l, max_len, rng)
    insertions += [pre_string[:l]]
    lengths += [l]
  return insertions, lengths


def parse_sequence(bytes seq, int k=10, kmers={}):
  """Go through the sequence filling out the k-mer dictionary

  :param seq:   the sequence
  :param k:     the k-mer length
  :param kmers: the kmer dictionary
  :return: (changes kmers in place)
  """
  cdef:
    char *c = seq
    int l = len(seq), n

  for n in xrange(l - k):
    try:
      kmers[seq[n:n+k]] += 1
    except KeyError:
      kmers[seq[n:n+k]] = 1


cdef float sequence_k_mer_score(bytes seq, int k, k_mer_count_table):
  cdef:
    int n, cnt
    float score = 0
    char *c = seq
  for n in range(len(seq) - k):
    score += k_mer_count_table.get(c[n:n + k], 1)
    cnt += 1
  return score / cnt


def score_long_sequence(bytes seq, int step, int k, k_mer_count_table):
  """Step through seq in step increments scoring each segment

  :param seq:
  :param step:
  :param k:
  :param k_mer_count_table:
  :return:
  """
  cdef:
    char *c = seq
    int i, idx = 0
  scores = np.empty(int(len(seq)/float(step)), dtype=np.uint32)
  for i in range(scores.shape[0]):
    scores[i] = k_mer_count_table.get(c[idx:idx + k], 1)
    idx += step
  return scores


def score_sequences_by_k_mer_count(sequence_list, k_mer_count_table):
  """Given a list of sequences and a k-mer count table score the sequence based on k-mer count content

  :param sequence_list:      list of sequences
  :param k_mer_count_table:  the k-mer count table dictionary. Keys = k-mers, values = counts in genome
  :return: list of scores, same size as sequence_list
  """
  if len(k_mer_count_table) == 0: return []
  k = len(k_mer_count_table.keys()[0])
  return [sequence_k_mer_score(seq, k, k_mer_count_table) for seq in sequence_list]
