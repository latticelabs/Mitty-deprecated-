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
  variant_waypoint = [(-1, -1, -1)]  # The start waypoint, guaranteed to be to the left and out of range of any base
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
        variant_waypoint += [(pos_ref, pos_alt, len(alt[c[0]]) - len(ref[c[0]]))]
        pos_alt += len(alt[c[0]])
      else:
        alt_fragments += [ref[c[0]]]
        pos_alt += len(ref[c[0]])
      pos_ref = stop[c[0]]
      c = next(c_iter, None)
  alt_fragments += [ref_seq[pos_ref:]]
  variant_waypoint += [(2 ** 31 - 1, 2 ** 31 - 1, 2 ** 31 - 1)]  # The end waypoint, guaranteed to be to the right of any base
  dtype = [('ref_pos', 'i4'), ('alt_pos', 'i4'), ('delta', 'i4')]
  return ''.join(alt_fragments), np.rec.fromrecords(variant_waypoint, dtype=dtype)


pref = lambda x, s: str(x) + s if x > 0 else ''

state_machine = {
  'b': {
    'S': lambda l, k, x: pref(l - 1, 'M') + '1X',
    'D': lambda l, k, x: pref(l, 'M'),
    'I': lambda l, k, x: pref(l, 'M'),
    'e': lambda l, k, x: pref(l, 'M')
  },
  'bI': {
    'S': lambda l, k, x: pref(k, 'S') + pref(l - 1 - k, 'M') + '1X',
    'D': lambda l, k, x: pref(k, 'S') + pref(l - k, 'M'),
    'I': lambda l, k, x: pref(k, 'S') + pref(l - k, 'M'),
    'e': lambda l, k, x: pref(k, 'S') + pref(l - k, 'M')
  },
  'S': {
    'S': lambda l, k, x: pref(l - 1, 'M') + '1X',
    'D': lambda l, k, x: pref(l, 'M'),
    'I': lambda l, k, x: pref(l, 'M'),
    'e': lambda l, k, x: pref(l, 'M')
  },
  'D' : {
    'S': lambda l, k, x: pref(x, 'D') + pref(l - 1, 'M') + '1X',
    'D': lambda l, k, x: pref(x, 'D') + pref(l, 'M'),
    'I': lambda l, k, x: pref(x, 'D') + pref(l, 'M'),
    'e': lambda l, k, x: pref(x, 'D') + pref(l, 'M')
  },
  'I': {
    'S': lambda l, k, x: pref(k, 'I') + pref(l - 1 - k, 'M') + '1X',
    'D': lambda l, k, x: pref(k, 'I') + pref(l - k, 'M'),
    'I': lambda l, k, x: pref(k, 'I') + pref(l - k, 'M'),
    'e': lambda l, k, x: pref(k, 'I') + pref(l - k, 'M'),
    'eI': lambda l, k, x: pref(l, 'S')
  }
}


def roll_cigars(variant_waypoints, reads):
  """Use beacons to generate POS and CIGAR strings for reads

  :param variant_waypoints: recarray, as returned by expand_sequence (pos_ref, pos_alt, delta)
  :param reads: numpy recarray with fields 'start_a' and 'read_len'
  :return: list of CIGAR strings same length as reads array
  """
  vw_r, vw_a, dl = variant_waypoints['ref_pos'], variant_waypoints['alt_pos'], variant_waypoints['delta']
  st_a, rd_len = reads['start_a'], reads['read_len']
  beacon_to_right = np.searchsorted(vw_a, st_a)
  cigars = []
  for n in range(reads.shape[0]):
    this_cigar = ''
    b1 = beacon_to_right[n]
    curr_pos = st_a[n]
    read_stop = st_a[n] + rd_len[n]
    state = 'b' if vw_a[b1 - 1] + dl[b1 - 1] < curr_pos else 'bI'
    while state != 'e' and state != 'eI':
      print state, curr_pos
      l = min(vw_a[b1] - curr_pos, rd_len[n])
      k = vw_a[b1 - 1] + dl[b1 - 1] - curr_pos
      x = -dl[b1 - 1]
      if vw_a[b1] > read_stop:
        new_state = 'eI' if state == 'I' else 'e'
      elif dl[b1] == 0:
        new_state = 'S'
      elif dl[b1] < 0:
        new_state = 'D'
      else:
        new_state = 'I'
      this_cigar += state_machine[state][new_state](l, k, x)
      curr_pos = vw_a[b1]
      state = new_state
      b1 += 1
    cigars += [this_cigar]
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