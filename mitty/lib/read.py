"""Declares a read structure that is commonly used."""
from ctypes import *


def _pp_seq(seq):
  return seq if len(seq) < 40 else seq[:20] + ' ... ' + seq[-20:]


class Read(Structure):
  _fields_ = [("POS", c_int32),
              ("CIGAR", c_char_p),
              ("perfect_seq", c_char_p),
              ("corrupted_seq", c_char_p),
              ("PHRED", c_char_p),  # Refers to the corrupted sequence
              ("_start_idx", c_int32),
              ("_stop_idx", c_int32)]  # internal use, refers to the pos_array. Used by roll_cigar etc

  def __repr__(self):
    return '(POS={0}, CIGAR="{1}", seq="{2}")'.format(self.POS, self.CIGAR, _pp_seq(self.perfect_seq))
