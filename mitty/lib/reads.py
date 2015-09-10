"""Infrastructure to handle reads"""
import re

import numpy as np


def expand_sequence(ref_seq, ml, chrom, copy):
  """Apply the variants in the list and return the consensus sequence

  :param ref_seq:    reference sequence
  :param ml:     master list of variants
  :param chrom:  [(no, het) ...] list of variants pointing to master list
                 no -> index on ml,
                 het -> 0 = copy 0, 1 = copy 1, 2 = homozygous
  :param copy:   0/1 which copy of the chromosome
  :return alt_seq, variant_waypoint, var_loc_alt_coordinates

  variant_waypoint -> recarray with the fields
      consensus sequence and array used by roll_cigars to determine POS and CIGAR strings for reads
          pos_ref: position on ref seq
          pos_alt: position on alt seq
          delta:  +k for insertions of length k, -k for deletions of length k, 0 for SNPs
  var_loc_alt_coordinates -> array of variants locations in the expanded sequence coordinates
  """
  pos_ref, pos_alt = 0, 0  # Current position in ref and alt coordinates
  alt_fragments = []
  variant_waypoint = [(-1, -1, 0)]  # The start waypoint, guaranteed to be to the left and out of range of any base and not an insertion or deletion
  var_loc_alt_coordinates = []
  pos, stop, ref, alt = ml.variants['pos'], ml.variants['stop'], ml.variants['ref'], ml.variants['alt']
  c_iter = chrom.__iter__()
  variant = next(c_iter, None)
  while variant is not None:
    if pos_ref < pos[variant[0]]:
      alt_fragments += [ref_seq[pos_ref:pos[variant[0]]]]
      pos_alt += pos[variant[0]] - pos_ref
      pos_ref = pos[variant[0]]
    else:
      if pos_ref == pos[variant[0]]:
        var_loc_alt_coordinates += [pos_alt]
        if variant[1] == 2 or variant[1] == copy:  # The variant applies to this chromosome copy
          alt_fragments += [alt[variant[0]]]
          dl = len(alt[variant[0]]) - len(ref[variant[0]])
          if dl == 0:
            variant_waypoint += [(pos_ref, pos_alt, dl)]  # For SNPs the waypoints don't move, so ref/alt stay same
          else:
            variant_waypoint += [(pos_ref + len(ref[variant[0]]), pos_alt + 1, dl)]
            # We shift the waypoint position to be the first non-match base
          pos_alt += len(alt[variant[0]])
        else:  # Skip this variant
          alt_fragments += [ref[variant[0]]]
          pos_alt += len(ref[variant[0]])
        pos_ref = stop[variant[0]]
        #pos_ref += len(ref[variant[0]])
      variant = next(c_iter, None)
  alt_fragments += [ref_seq[pos_ref:]]

  final_delta = variant_waypoint[-1][0] - variant_waypoint[-1][1]
  if final_delta > 0:
    final_ref, final_alt = 2 ** 31 - 1, 2 ** 31 - 1 - final_delta
  else:
    final_ref, final_alt = 2 ** 31 - 1 - final_delta, 2 ** 31 - 1
  variant_waypoint += [(final_ref, final_alt, -1)]
  # The end waypoint, guaranteed to be to the right of any base and not a SNP, and maintaining the delta
  dtype = [('ref_pos', 'i4'), ('alt_pos', 'i4'), ('delta', 'i4')]
  return ''.join(alt_fragments), np.rec.fromrecords(variant_waypoint, dtype=dtype), var_loc_alt_coordinates


# TODO: make this code more elegant
def roll_cigars(variant_waypoints, reads):
  """Use beacons to generate POS and CIGAR strings for reads

  :param variant_waypoints: recarray, as returned by expand_sequence (pos_ref, pos_alt, delta)
  :param reads: numpy recarray with fields 'start_a' and 'read_len'
  :return: pos, cigars
     - list of POS values
     - list of CIGAR strings same length as reads array
  """
  v_r, v_a, dl = variant_waypoints['ref_pos'], variant_waypoints['alt_pos'], variant_waypoints['delta']
  rd_st, rd_len = reads['start_a'], reads['read_len']
  waypoint_right = np.searchsorted(v_a, rd_st)
  cigars = []
  pos = []
  for rd_no in range(reads.shape[0]):
    r_start = rd_st[rd_no]
    r_stop = rd_st[rd_no] + rd_len[rd_no] - 1
    n = waypoint_right[rd_no]

    m = min(v_a[n], r_stop + 1) - r_start
    cigar = str(m) + '=' if m > 0 else ''
    this_pos = v_r[n - 1] + r_start - v_a[n - 1]  # In our system the previous waypoint has the delta between ref and alt
    if dl[n - 1] > 0:  # The previous variant was an insertion, possibility for soft-clipping
      this_pos = v_r[n - 1] + max(r_start - v_a[n - 1] - dl[n - 1], 0)  # POS is nearest match base (and INS has been shifted one base to the right for waypoint)
      sc = v_a[n - 1] + dl[n - 1] - r_start
      if sc > 0:  # Yes, a soft-clip
        cigar = str(sc) + 'S' + (str(m - sc) + '=' if m - sc > 0 else '')
    if r_start == v_a[n] and dl[n] < 0:  # Corner case: we are starting at a deletion
      this_pos = v_r[n] + r_start - v_a[n]

    pos += [this_pos]  # POS is 0 indexed as per BAM spec
    while r_stop >= v_a[n]:
      if dl[n] == 0:  # SNP
        m = min(v_a[n+1], r_stop + 1) - v_a[n] - 1
        cigar += '1X' + (str(m) + '=' if m > 0 else '')
      elif dl[n] > 0:  # INS
        if r_start == v_a[n]:  # Corner case: we are starting right at an insertion
          if v_a[n] + dl[n] - 1 >= r_stop:  # Completely inside insertion
            cigar = str(rd_len[rd_no]) + 'S'
          else:  # Soft-clipped, then with Ms
            m = min(v_a[n + 1], r_stop + 1) - r_start
            sc = min(dl[n], r_stop + 1 - r_start)
            cigar = str(sc) + 'S' + str(m - sc) + '='
        elif v_a[n] + dl[n] - 1 < r_stop:  # Insert has anchor on other side
          m = min(v_a[n + 1], r_stop + 1) - v_a[n] - dl[n]
          cigar += str(dl[n]) + 'I' + (str(m) + '=' if m > 0 else '')
        else:  # Handle soft-clip at end
          cigar += str(r_stop - v_a[n] + 1) + 'S'
      else:  # DEL
        m = min(v_a[n + 1], r_stop + 1) - v_a[n]
        if r_start != v_a[n]:
          cigar += str(-dl[n]) + 'D' + str(m) + '='
        else:
          cigar += str(m) + '='  # Corner case: if we start right at a deletion
      n += 1
    cigars += [cigar]
  return pos, cigars


cig_re = re.compile(r'(\d+?)M(\d+?)M')


def old_style_cigar(cigar):
  """Given an extended cigar ('X's for mismatch and '=' for match) convert it into an old style cigar with 'M's and no
  'X's

  :param cigar:
  :return: old style cigar
  """
  cigar = cigar.replace('=', 'M')  # First replace all '='s with 'M's
  if 'X' in cigar:  # More complicated - need to merge Xs and =s into Ms as needed
    cigar = cigar.replace('X', 'M')
    # Now we need to collapse the 'M's as needed
    cigar, n = cig_re.subn(lambda m: str(int(m.group(1)) + int(m.group(2))) + 'M', cigar)
    while n:
      cigar, n = cig_re.subn(lambda m: str(int(m.group(1)) + int(m.group(2))) + 'M', cigar)
  return cigar
