"""Infrastructure to handle reads"""
import numpy as np


def expand_sequence(ref_seq, ml, chrom, copy):
  """Apply the variants in the list and return the consensus sequence

  :param ref_seq:    reference sequence
  :param ml:     master list of variants
  :param chrom:  list of variants, referring to master list
  :param copy:   0/1 which copy of the chromosome
  :return alt_seq, beacons: consensus sequence and array used by roll_cigars to determine POS and CIGAR strings for reads

  beacons -> recarray with the fields
              pos_ref: position on ref seq
              pos_alt: position on alt seq
              delta:  +k for insertions of length k, -k for deletions of length k, 0 for SNPs
  """
  pos_ref, pos_alt = 0, 0  # Current position in ref and alt coordinates
  alt_fragments = []
  beacons = []
  pos, stop, ref, alt = ml.variants['pos'] - 1, ml.variants['stop'] - 1, ml.variants['ref'], ml.variants['alt']
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
        beacons += [(pos_ref, pos_alt, len(alt[c[0]]) - len(ref[c[0]]))]
        pos_alt += len(alt[c[0]])
      else:
        alt_fragments += [ref[c[0]]]
        pos_alt += len(ref[c[0]])
      pos_ref = stop[c[0]]
      c = next(c_iter, None)
  alt_fragments += [ref_seq[pos_ref:]]
  dtype = [('ref_pos', 'i4'), ('alt_pos', 'i4'), ('delta', 'i4')]
  return ''.join(alt_fragments), np.rec.fromrecords(beacons, dtype=dtype) if len(beacons) else np.rec.recarray((0, 3), dtype=dtype)


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