"""This module contains functions that generate simulated reads. The module can be called as a script as well.

Commandline::

  Usage:
    vcf2reads [localreads]  --wg=WG  --vcf=VCF  --out=OUT  --paramfile=PFILE  [--corrupt]  [--fastq] [--read_border_size=RBS] [--reads_per_block=BL] [--master_seed=MS] [-v]

  Options:
    localreads              Take reads around insertions only
    --wg=WG                 Whole genome file
    --vcf=VCF               VCF file
    --out=OUT               Output file name prefix
    --read_border_size=RBS  How many bases before and after the insertion do we include in our window [default: 500]
    --paramfile=PFILE       Name for parameter file
    --corrupt               Write out corrupted reads too.
    --fastq                 Write as FASTQ instead of BAM (simulated_reads.fastq)
    --reads_per_block=BL    Generate these many reads at a time (Adjust to machine resources). [default: 100000]
    --master_seed=MS        If this is specified this passes a master seed to the read plugin.
                            This overrides any individual seeds specified by the parameter file.
    -v                      Dump detailed logger messages

Parameter file example::

  {
    "take reads from": [1,2],
    "coverage": 5.0,
    "output_file_prefix": "sim_reads",
    "read_model": "simple_reads",
    "model_params": {
      "paired": false,
      "read_len": 100,
      "template_len": 250
    }
  }

Note:
Expanding the variant sequence takes a bunch of memory because we need 14 bytes per base:

1 byte - forward seq
1 byte - complement seq
4 bytes - pos array  (needed for POS)
4 bytes - diff pos array (needed for CIGAR)
4 bytes - offset array (needed for offset for reads deep in inserts)
-
14 bytes

For this reason we adopt a just in time expansion where by we feed chunks of the variant sequence to the read plugin.
The smaller the chunk size the less extra memory is needed but the chunk computation has an overhead. In general, you
should make the chunk size as large as you can given your memory


Algorithm:

1. Read plugin gives us a list of read + template positions (sorted by position along the mutated sequence).
2. We start at the beginning of the reference sequence
3. The read reference is the reference sequence
4. If the template end is before the next variant position:
   1. take reads, repeat 4.
5. If the template end crosses the next variant position:
   1. If read reference is the reference sequence, start an alt sequence,
      otherwise clip existing alt sequence (and pos_array)
   2. Expand the variant, add it to alt seq
   3. Repeat 2 as needed to go past template
   4. Goto 4
6. Repeat all this until reads are exhausted


Roadmap

1. Rewrite to use only reference and VCF file to generate reads
2. Think about how to separate out plugin from main code etc.
3. Make corruption plugin, read plugin separate?


1. The quality scores are in Phred scale (as specified in the SAM spec)
2. We supply the prefix of output file name in the parameter file . Say we set this as sim_reads.
   The perfect reads will be saved to sim_reads.bam (or sim_reads.fastq). If we ask for corrupted reads
   we will get the corrupted reads in the file sim_reads_c.fastq.
   A text sidecar file sim_reads.info will always be saved with simulation parameters.

3 *** Can probably refactor the whole file for better readability **
  *** see if we can write shorter functions, see if we can reduce the number of parameters/bundle them ***
  *** further efficiency gains will probably be minimal - the bottle neck function has been cythonized ***

"""
import numpy
from mitty.variation import *  # Yes, it's THAT important
import string

DNA_complement = string.maketrans('ATCGN', 'TAGCN')
pos_null = numpy.empty((0,), dtype='u4')  # Convenient, used in apply_one_variant


class Read(Structure):
  _fields_ = [("POS", c_int32),
              ("CIGAR", c_char_p),
              ("perfect_seq", c_char_p),
              ("corrupted_seq", c_char_p),
              ("PHRED", c_char_p)]  # Refers to the corrupted sequence

  def __eq__(self, other):
    # This does not do an isinstance check for speed reasons.
    if self.POS == other.POS and self.CIGAR == other.CIGAR and self.seq == other.seq and self.PHRED == other.PHRED:
      return True
    else:
      return False

  def __ne__(self, other):
      return not self.__eq__(other)

  def __repr__(self):
    return '(POS={0}, CIGAR="{1}", seq="{2}")'.format(self.POS, self.CIGAR,
                                                      self.seq if len(self.seq) < 40 else self.seq[:20] +
                                                      ' ... ' + self.seq[-20:])


def init_read_model(read_model):
  return None


def get_variant_sequence_generator(ref_chrom_seq='', c1=[], chrom_copy=0, block_len=10e6, over_lap_len=200):
  """Return the computed variant sequence in blocks
  Args:
    ref_chrom_seq   : (string) the reference sequence
    c1              : (list of Variant) standard chromosome format
    chrom_copy      : Which copy of the chromosome (0 or 1)
    block_len       : how many bases of the variant sequence to return per iteration
    overlap_len     : how many bases of the old sequence to include in the new block
                      This should be hinted by the expected template length of the read model
  Returns:
    (start_idx      : (u4) in variant sequence coordinates - where does this block start
     seq, c_seq,    : (strings) variant sequence and complement
     arr            : (list of numpy.array u4) [a1, a2, a3  ]
                      a1 -  position array used for POS
                      a2 -  diff(p_arr) used for CIGARs
                      a3 -  offset array used for offset calculation for reads deep in inserts
  Raises:
    StopIteration   : When we are all done

  Notes:
    1. This function has not been refactored to avoid function overhead
    2. We use a list of arrays for arr rather than a 2D array as 2D slicing turns out to be more expensive
  """
  l_ref_seq = len(ref_chrom_seq)
  cc = chrom_copy

  over_lap = 0
  seq_fragments = ['']
  arr_null = numpy.empty(0, dtype='u4')  # Convenient
  arr_fragments = [[arr_null for _ in [0, 1, 2]]]
  ptr, var_ptr = 0, 0
  var_ptr_start = var_ptr
  var_ptr_finish = var_ptr + block_len
  c1_iter = c1.__iter__()
  variant = next(c1_iter, None)
  while ptr < l_ref_seq:
    if variant is None or (ptr < variant.POS - 1):  # We should copy just the reference
      if variant is None: # No more variants left
        ref_ptr_start = ptr
        ptr = min(ptr + var_ptr_finish - var_ptr, l_ref_seq)
        var_ptr += ptr - ref_ptr_start
      else:  # Copy as much of the reference as we can
        ref_ptr_start = ptr
        ptr = min(ptr + var_ptr_finish - var_ptr, variant.POS - 1)
        var_ptr += ptr - ref_ptr_start
      seq_fragments += [ref_chrom_seq[ref_ptr_start:ptr]]
      arr_fragments += [[
        numpy.arange(ref_ptr_start + 1, ptr + 1, dtype='u4'),
        numpy.empty(ptr - ref_ptr_start, dtype='u4'),
        numpy.zeros(ptr - ref_ptr_start, dtype='u4')
      ]]
    else:  # variant exists and we are on it, expand it
      alt, ref = variant.ALT, variant.REF
      if (variant.het == HET1 and cc == 1) or (variant.het == HET2 and cc == 0):
        alt = variant.REF

      l_alt, l_ref = len(alt), len(ref)
      ptr_adv = l_ref  # When we get to this function, our pointer is sitting at POS
      pos_alt = numpy.arange(variant.POS, variant.POS + min(l_alt, l_ref), dtype='u4')  # We might have ref bases in alt
      if l_alt > l_ref:  # This was an insertion
        pos_alt = numpy.concatenate((pos_alt, numpy.ones(l_alt - l_ref, dtype='u4') * variant.stop))

      seq_fragments += [alt]
      arr_fragments += [[pos_alt, numpy.empty(l_alt, dtype='u4'), numpy.arange(l_alt, dtype='u4')]]

      ptr += ptr_adv
      var_ptr += l_alt
      variant = next(c1_iter, None)  # Load the next variant in preparation

    if var_ptr >= var_ptr_finish or ptr >= l_ref_seq:  # Ok, we've got enough for this block
      # Data to return
      this_idx = var_ptr_start - over_lap
      this_seq_block = ''.join(seq_fragments)
      this_c_seq_block = this_seq_block.translate(DNA_complement)
      this_arr = [numpy.concatenate([a_frag[n] for a_frag in arr_fragments]) for n in [0, 1, 2]]
      this_arr[1][:-1] = numpy.diff(this_arr[0])
      this_arr[1][-1] = (ptr + 1) - this_arr[0][-1]

      # Prep for next block before yield, there is no coming back
      var_ptr_start = var_ptr
      var_ptr_finish = var_ptr + block_len
      over_lap = min(over_lap_len, len(this_seq_block))
      seq_fragments = [this_seq_block[-over_lap:]]
      arr_fragments = [[a_frag[-over_lap:] for a_frag in this_arr]]

      yield this_idx, this_seq_block, this_c_seq_block, this_arr


def roll_cigar(pos_array, pos=0):
  """Given a position array write out a POS and CIGAR value."""
  if numpy.count_nonzero(pos_array[1]) == 0:  #This is a read from the middle of an insertion
    align_pos = pos_array[0][0] - 1  # This assumes we use the A -> ATCG form for VCF lines
    cigar = '>' + str(pos_array[2][0])  # The offset


  mapped = False
  cigar = ''
  counter = 0
  cigar_fragment = None
  for n in range(pos_array.size):
    dp = pos_array[n+1] - pos_array[n]
    if dp == 1:
      mapped = True  # As long as we have one M we are a mapped read
      if cigar_fragment != 'M':
        if counter > 0:  # Flush
          cigar += '{:d}{:s}'.format(counter, cigar_fragment)
          counter = 0
      cigar_fragment = 'M'
      counter += 1
    elif dp == 0:
      if cigar_fragment != 'I':
        if counter > 0:  # Flush
          cigar += '{:d}{:s}'.format(counter, cigar_fragment)
          counter = 0
      cigar_fragment = 'I'
      counter += 1
    elif dp > 1:
      mapped = True  # As long as we have one M we are a mapped read
      if cigar_fragment != 'M':
        if counter > 0:  # Flush
          cigar += '{:d}{:s}'.format(counter, cigar_fragment)
          counter = 0
      cigar_fragment = 'M'  # We need to set this because we could be at the start of a read and type = None still
      counter += 1
      cigar += '{:d}{:s}'.format(counter, cigar_fragment)
      cigar_fragment = 'D'
      counter = dp - 1

  if cigar_fragment != 'D':  # Flush all but 'D'. We only write D if we cross a D boundary
    cigar += '{:d}{:s}'.format(counter, cigar_fragment)

  if mapped:
    align_pos = pos_array[0]
  else:
    align_pos = pos  # Unmapped read, deep in the middle of an insert. This is where we are!
    cigar = '>' + str(offset)  # A deep CIGAR, indicates last POS (has to be insert) and offset into insert

  return align_pos, cigar


# def perfect_reads_from_seq(ref=[], pos_array=None, offset_pos_array=None, offset_array=None, read_info=[], counts=[]):
#   """
#   Args:
#     read_info      : list of [(pos, len, dir, strand) ... ]
#                       e.g [(pos, len, dir, strand), (pos2, len2, dir2, stand2)] if paired
#                      pos = start of read
#                      len = len of read
#                      dir = +1 for forward, -1 for reverse
#     counts         : How many of each template do we have
#   """
#   # phred = '~' * max_read_len
#   # perfect_cigar = '{:d}M'.format(max_read_len)
#   read_list = []
#   #append = read_list.append
#   for this_template_info, repeat_count in zip(read_info, counts):
#     if repeat_count == 0: continue  # No reads here
#     for (r_pos, r_len, r_dir, r_strand) in this_template_info:
#       if pos_array is not None: # This sequence contains variants
#         pos, cigar = roll_cigar(pos_array[r_pos:r_pos + r_len], offset=offset_array[x1])
#       else:  # This is part of the reference sequence
#         pos1, cigar1 = x1, perfect_cigar  # POS is 1 indexed
#
#     if paired:
#       x1, x2 = ref_start + n + template_len - 1, ref_start + n + template_len - read_len - 1
#       if pos_array is not None:
#         pos2, cigar2 = roll_cigar(pos_array[x1:x2:-1], offset=offset_array[x2])
#       else:
#         pos2, cigar2 = x2 + 1, perfect_cigar
#
#     for strand in strands:
#       this_read_seq = ref[strand][ref_start + n:ref_start + n + read_len]
#       this_read = [[Read(pos1, cigar1, this_read_seq, phred)]]
#       if paired:
#         that_read_seq = ref[1 - strand][ref_start + n + template_len - 1:ref_start + n + template_len - read_len - 1:-1]
#         this_read = [[this_read[0][0], Read(pos2, cigar2, that_read_seq, phred)]]
#       read_list += this_read
#   return read_list


# def expand_variant_seq():
#
#
#
#
# def perfect_reads_from_chrom_copy(ref_chrom_seq='', c1=[],
#                                   template_len_hint=None,
#                                   block_len_hint=None,
#                                   read_pos_iter=None,
#                                   read_iter=None):
#   read_list = []
#   ref_pointer = 0  # Where are we on the reference sequence
#   alt_pointer = 0  # Where are we on the variant sequence
#   c1_iter = c1.__iter__()
#   alt_chrom_seq = ''
#   variant = next(c1_iter, None)
#   next_template_pos = next(read_pos_iter, None)
#   while 1:
#     while next_template_pos
#
#
#
#   for variant in c1:
#     # Take reads up to this variant
#     if ref_pointer + template_len_hint < variant.POS: #
#
#
#
#
#
#     copy[n].append(ref_seq[pointer[n]:variant.POS - 1])
#     pos[n].append(numpy.arange(pointer[n] + 1, variant.POS, dtype='u4'))
#     pointer[n] += max(0, variant.POS - 1 - pointer[n])
#     mutated_pointer[n] += len(copy[n][-1])
#
#
#     read_list += perfect_reads_from_ref_seq(chrom_seq[ref_pointer:variant.POS])
#
#
#
# def perfect_reads_from_genome(ref={}, g1={}, chrom_list=[], p=0.1, k=100):
#   """
#   Args:
#     ref        : reference genome
#     g1         : genome as retrieved from vcf file
#     chrom_list : list of chromosomes to take reads from
#     p          : probability of read per base (bernoulli parameters for reads)
#     k          : number of 'tries' per base
#
#   Returns:
#     read_list  : a list of Read values
#   """
#   read_list = []
#   for chrom in chrom_list:
#     read_list += perfect_reads_from_chrom(ref[chrom], g1[chrom], p=p, k=k)
#   return read_list
#
#
# def corrupt_reads(read_list, corruption_model={}):
#   pass

def generate_reads(ref={}, g1={}, chrom_list=[], read_model=None, block_size=100):
  """This is written as a generator because we might have a lot of reads and we want to flush the reads to file as
  we go.

  Args:
    ref        : reference genome
    g1         : genome as retrieved from vcf file
    chrom_list : list of chromosomes to take reads from
    read_model :
    block_size : How many bases do we send to the read generator

  Returns:
    read_gen   : A generator that will keep giving reads until we are done

  1. Read plugin gives us a list of read + template positions (sorted by position along the mutated sequence).
  2. We start at the beginning of the reference sequence
  3. The read reference is the reference sequence
  4. If the template end is before the next variant position:
     1. take reads, repeat 4.
  5. If the template end crosses the next variant position:
     1. If read reference is the reference sequence, start an alt sequence,
        otherwise clip existing alt sequence (and pos_array)
     2. Expand the variant, add it to alt seq
     3. Repeat 2 as needed to go past template
     4. Goto 4
  6. Repeat all this until reads are exhausted

  """
  read_iter = read_model.initialize(block_size=100)
  for read in read_iter:
    pass





