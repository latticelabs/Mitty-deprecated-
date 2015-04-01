"""Infrastructure to handle reads"""
import numpy as np


def expand_sequence(ref_seq, ml, chrom, copy):
  """Apply the variants in the list and return the consensus sequence

  :param ref_seq:    reference sequence
  :param ml:     master list of variants
  :param chrom:  list of variants, referring to master list
  :param copy:   0/1 which copy of the chromosome
  :return alt_seq, variant_waypoint: consensus sequence and array used by roll_cigars to determine POS and CIGAR strings for reads

  variant_waypoint -> recarray with the fields
              pos_ref: position on ref seq
              pos_alt: position on alt seq
              delta:  +k for insertions of length k, -k for deletions of length k, 0 for SNPs

  """
  pos_ref, pos_alt = 0, 0  # Current position in ref and alt coordinates
  alt_fragments = []
  variant_waypoint = [(-1, -1, 0)]  # The start waypoint, guaranteed to be to the left and out of range of any base and not an insertion or deletion
  pos, stop, ref, alt = ml.variants['pos'], ml.variants['stop'], ml.variants['ref'], ml.variants['alt']
  c_iter = chrom.__iter__()
  c = next(c_iter, None)
  while c:
    if pos_ref < pos[c[0]]:
      alt_fragments += [ref_seq[pos_ref:pos[c[0]]]]
      pos_alt += pos[c[0]] - pos_ref
      pos_ref = pos[c[0]]
    else:
      if c[1] == 2 or c[1] == copy:
        alt_fragments += [alt[c[0]]]
        dl = len(alt[c[0]]) - len(ref[c[0]])
        dp = 0 if dl == 0 else 1  # We shift the waypoint position to be the first non-match base
        variant_waypoint += [(pos_ref + dp, pos_alt + dp, dl)]
        pos_alt += len(alt[c[0]])
      else:
        alt_fragments += [ref[c[0]]]
        pos_alt += len(ref[c[0]])
      pos_ref = stop[c[0]]
      c = next(c_iter, None)
  alt_fragments += [ref_seq[pos_ref:]]
  variant_waypoint += [(2 ** 31 - 1, 2 ** 31 - 1, -1)]  # The end waypoint, guaranteed to be to the right of any base and not a SNP
  dtype = [('ref_pos', 'i4'), ('alt_pos', 'i4'), ('delta', 'i4')]
  return ''.join(alt_fragments), np.rec.fromrecords(variant_waypoint, dtype=dtype)


# TODO: make this code more elegant
def roll_cigars(variant_waypoints, reads):
  """Use beacons to generate POS and CIGAR strings for reads

  :param variant_waypoints: recarray, as returned by expand_sequence (pos_ref, pos_alt, delta)
  :param reads: numpy recarray with fields 'start_a' and 'read_len'
  :return: list of CIGAR strings same length as reads array
  """
  v_r, v_a, dl = variant_waypoints['ref_pos'], variant_waypoints['alt_pos'], variant_waypoints['delta']
  rd_st, rd_len = reads['start_a'], reads['read_len']
  waypoint_right = np.searchsorted(v_a, rd_st)
  cigars = []
  for rd_no in range(reads.shape[0]):
    r_start = rd_st[rd_no]
    r_stop = rd_st[rd_no] + rd_len[rd_no] - 1
    n = waypoint_right[rd_no]

    m = min(v_a[n], r_stop + 1) - r_start
    cigar = str(m) + 'M' if m > 0 else ''
    if dl[n - 1] > 0:  # The previous variant was an insertion, possibility for soft-clipping
      sc = v_a[n - 1] + dl[n - 1] - r_start
      if sc > 0:  # Yes, a soft-clip
        cigar = str(sc) + 'S' + (str(m - sc) + 'M' if m - sc > 0 else '')
    if dl[n] > 0 and r_start == v_a[n]:  # Corner case: we are starting right at an insertion
      if v_a[n] + dl[n] - 1 >= r_stop:  # Completely inside insertion
        cigar = str(rd_len[rd_no]) + 'S'
      else:  # Soft-clipped, then with Ms
        m = min(v_a[n + 1], r_stop + 1) - r_start
        sc = min(dl[n], r_stop + 1 - r_start)
        cigar = str(sc) + 'S' + str(m - sc) + 'M'
      n += 1
    while r_stop >= v_a[n]:
      if dl[n] == 0:  # SNP
        m = min(v_a[n+1], r_stop + 1) - v_a[n] - 1
        cigar += '1X' + (str(m) + 'M' if m > 0 else '')
      elif dl[n] > 0:  # INS
        if v_a[n] + dl[n] - 1 < r_stop:  # Insert has anchor on other side
          m = min(v_a[n + 1], r_stop + 1) - v_a[n] - dl[n]
          cigar += str(dl[n]) + 'I' + (str(m) + 'M' if m > 0 else '')
        else:  # Handle soft-clip at end
          cigar += str(r_stop - v_a[n] + 1) + 'S'
      else:  # DEL
        cigar += str(-dl[n]) + 'D' + str(v_a[n+1] - v_a[n]) + 'M'
      n += 1
    cigars += [cigar]
  return cigars


from ctypes import *
direction = '><'


def _pp_seq(seq):
  return seq if len(seq) < 40 else seq[:20] + ' ... ' + seq[-20:]


class Read(Structure):
  _fields_ = [("direction", c_char),      # Is the read forward '>' or reverse '<'
              ("POS", c_int32),
              ("CIGAR", c_char_p),
              ("perfect_seq", c_char_p),
              ("corrupt_seq", c_char_p),
              ("PHRED", c_char_p),        # Refers to the corrupted sequence
              ("_start_idx", c_int32),
              ("_stop_idx", c_int32)]     # internal use, refers to the pos_array. Used by roll_cigar etc

  def __repr__(self):
    return '(POS={0}, CIGAR="{1}", seq="{2}")'.format(self.POS, self.CIGAR, _pp_seq(self.perfect_seq))