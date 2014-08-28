"""This module contains functions that generate simulated reads. The module can be called as a script as well.

Commandline::

  Usage:
    vcf2reads  --fa_dir=FADIR  [--vcf=VCF]  --out=OUT  --pfile=PFILE  [--corrupt]  [--block_len=BL] [--master_seed=MS] [-v]

  Options:
    --fa_dir=FADIR          Directory where genome is located
    --vcf=VCF               VCF file. If not given we will take reads from the reference genome
    --out=OUT               Output file name prefix
    --pfile=PFILE           Name for parameter file
    --corrupt               Write out corrupted reads too.
    --block_len=BL          Consider the sequence in chunks this big. See notes [default: 1000000]
    --master_seed=MS        If this is specified this passes a master seed to the read plugin.
                            This overrides any individual seeds specified by the parameter file.
                            If this is not set the random number generator will be initialized via the default method
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
should make the chunk size as large as you can given your memory keeping in mind that each base takes 14 bytes of space
when the variant sequence is expanded.


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
import docopt
import numpy
import vcf
import json
from mitty.variation import *  # Yes, it's THAT important
from mitty.utility.genome import FastaGenome
import string
DNA_complement = string.maketrans('ATCGN', 'TAGCN')
pos_null = numpy.empty((0,), dtype='u4')  # Convenient, used in apply_one_variant
int2str = [str(n) for n in range(1001)]  # int2str[n] is faster than str(n).
# the read model should give max read length and we should use init_int2str to (re)generate an appropriately sized table
# Leave this initialization in as it makes some testing easier


def init_int2str(max_read_len):
  global int2str
  if max_read_len >= len(int2str):
    int2str = [str(n) for n in range(max_read_len + 1)]


import logging
logger = logging.getLogger(__name__)

SEED_MAX = 10000000000  # For each call of the JIT expander we pass a random seed.
                        # We generate these seeds from the given master seed and


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

  # def __eq__(self, other):
  #   # This does not do an isinstance check for speed reasons.
  #   if self.POS == other.POS and self.CIGAR == other.CIGAR and self.perfect_seq == other.perfect_seq and self.PHRED == other.PHRED:
  #     return True
  #   else:
  #     return False
  #
  # def __ne__(self, other):
  #     return not self.__eq__(other)

  def __repr__(self):
    return '(POS={0}, CIGAR="{1}", seq="{2}")'.format(self.POS, self.CIGAR, _pp_seq(self.perfect_seq))


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
    iter            : an iterator that yields tuples

     start_idx      : (u4) in variant sequence coordinates - where does this block start
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


#TODO handle softclip at start of read
def roll_cigar(pos_array, start_idx, stop_idx):
  cigar = ''
  counter = 0
  cigar_fragment = None
  for dp in pos_array[1][start_idx:stop_idx]:
    if dp == 1:
      if cigar_fragment != 'M':
        if counter > 0:  # Flush
          cigar += int2str[counter] + cigar_fragment  # '{:d}{:s}'.format(counter, cigar_fragment)
          counter = 0
      cigar_fragment = 'M'
      counter += 1
    elif dp == 0:
      if cigar_fragment != 'I':
        if counter > 0:  # Flush
          cigar += int2str[counter] + cigar_fragment  # '{:d}{:s}'.format(counter, cigar_fragment)
          counter = 0
      cigar_fragment = 'I'
      counter += 1
    elif dp > 1:
      if cigar_fragment != 'M':
        if counter > 0:  # Flush
          cigar += int2str[counter] + cigar_fragment  # '{:d}{:s}'.format(counter, cigar_fragment)
          counter = 0
      cigar_fragment = 'M'  # We need to set this because we could be at the start of a read and type = None still
      counter += 1
      cigar += int2str[counter] + cigar_fragment  # '{:d}{:s}'.format(counter, cigar_fragment)
      cigar_fragment = 'D'
      counter = dp - 1

  # Flush all but 'D'. We only write D if we cross a D boundary
  if cigar_fragment == 'M':
    cigar += int2str[counter] + cigar_fragment  # '{:d}{:s}'.format(counter, cigar_fragment)
  elif cigar_fragment == 'I':  # Ooh sophisticated - we write softclips now
    cigar += int2str[counter] + 'S'

  return cigar


def package_reads(template_list, pos_array):
  """Fills out POS and CIGAR in place"""
  for template in template_list:
    for read in template:
      start_idx, stop_idx = read._start_idx, read._stop_idx
      if numpy.any(pos_array[1][start_idx:stop_idx] != 1):  # Somewhere we have an indel
        if numpy.count_nonzero(pos_array[1][start_idx:stop_idx]) == 0:  #This is a read from the middle of an insertion
          align_pos = pos_array[0][start_idx] - 1  # This assumes we use the A -> ATCG form for VCF lines
          cigar = '>' + str(pos_array[2][start_idx])  # The offset
        else:  # Tough luck, gotta compute the CIGAR manually
          align_pos, cigar = pos_array[0][start_idx], roll_cigar(pos_array, start_idx, stop_idx)
      else:  # This is a perfect match
        align_pos, cigar = pos_array[0][start_idx], str(stop_idx - start_idx) + 'M'
      read.POS, read.CIGAR = align_pos, cigar


def reads_from_genome(ref={}, g1={}, chrom_list=[], read_model=None, model_params={}, block_len=10e6, master_seed=1):
  """
  Args:
    ref          : reference genome
    g1           : genome as retrieved from vcf file
    chrom_list   : list of chromosomes to take reads from
    read_model   : read model
    model_params :
    block_len    : expand block length
    master_seed  :

  Returns:
    iter         : a generator that returns

    read_list    : a list of completed Read objects. Paired reads come sequentially
    chrom        : chromosome the read was taken from [1,2....]
    cc           : chromosome copy the read was taken from [0,1]

  Raises:
    StopIteration   : When we are all done

  Notes:
    This is written as a generator because we might have a lot of reads and we want to flush the reads to file as we go.
  """
  seed_rng = numpy.random.RandomState(seed=master_seed)

  read_model_state = read_model.initialize(model_params, master_seed)
  # This is meant for us to store any precomputed tables/constants. If we wish we can store RNGs here
  overlap_len = read_model.overlap_len(read_model_state)
  max_read_len = read_model.max_read_len(read_model_state)
  # This is used to determine the overlap for the blocks we feed to the read generator
  init_int2str(max_read_len)  # And to compute this table, of course

  for chrom in chrom_list:
    if chrom not in ref: continue
    for cc in [0, 1]:
      vsg = get_variant_sequence_generator(ref_chrom_seq=ref[chrom], c1=g1.get(chrom, []), chrom_copy=cc,
                                           block_len=block_len, over_lap_len=overlap_len)
      for this_idx, this_seq_block, this_c_seq_block, this_arr in vsg:
        tl, read_model_state = read_model.generate_reads(this_idx, this_seq_block, this_c_seq_block, this_arr,
                                                         read_model_state, seed_rng.randint(SEED_MAX))
        # Note that we reseed the generator each time. Each chunk is assumed to be independent of the last
        # The model could ignore this independence by storing its RNGs in read_model_state
        package_reads(tl, this_arr)
        yield tl, chrom, cc


def write_reads_to_file(fp, fp_c, template_list, chrom, cc, serial_no, write_corrupted=False):
  """
  Args:
    fp           : file pointer
    fp_c         : file pointer for corrupted reads
    read_list    : a list of completed Read objects. Paired reads come sequentially
    chrom        : chromosome the read was taken from [1,2....]
    cc           : chromosome copy the read was taken from [0,1]
    serial_no    : serial number of first read
    corrupt      : write corrupted reads
  """
  # qname format is chrom:copy|rN|POS1|CIGAR1|POS2|CIGAR2
  for n, template in enumerate(template_list):
    paired = len(template) == 2
    qname = chrom.__str__() + ':' + cc.__str__() + '|r' + (n + serial_no).__str__() + '|' + template[0].POS.__str__() \
      + '|' + template[0].CIGAR
    if paired:
      qname += '|' + template[1].POS.__str__() + '|' + template[1].CIGAR

    fp.write('@' + qname + '\n' + template[0].perfect_seq + '\n+\n' + '~' * len(template[0].perfect_seq) + '\n')
    if paired:
      fp.write('@' + qname + '\n' + template[1].perfect_seq + '\n+\n' + '~' * len(template[1].perfect_seq) + '\n')

    if write_corrupted:
      fp_c.write('@' + qname + '\n' + template[0].corrupted_seq + '\n+\n' + template[0].PHRED + '\n')
      if paired:
        fp_c.write('@' + qname + '\n' + template[1].corrupted_seq + '\n+\n' + template[1].PHRED + '\n')


def main(fp, fp_c=None, ref={}, g1={}, chrom_list=[], read_model=None, model_params={}, block_len=10e6, master_seed=1):
  read_gen = reads_from_genome(ref=ref, g1=g1, chrom_list=chrom_list,
                               read_model=read_model, model_params=model_params,
                               block_len=block_len, master_seed=master_seed)
  write_corrupted = False if fp_c is None else True
  serial_no = 0
  for template_list, chrom, cc in read_gen:
    write_reads_to_file(fp, fp_c, template_list=template_list, chrom=chrom, cc=cc,
                        serial_no=serial_no, write_corrupted=write_corrupted)
    serial_no += len(template_list)


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  if args['-vv']:
    logging.basicConfig(level=logging.DEBUG)
  else:
    logging.basicConfig(level=logging.WARNING)

  if args['-v']:
    logger.setLevel(logging.DEBUG)

  fp = open(args['--out'] + '.fq')
  fp_c = open(args['--out'] + '_c.fq') if args['--corrupt'] else None
  ref_genome = FastaGenome(seq_dir=cmd_args['--fa_dir'])
  params = json.load(open(args['--pfile'], 'r'))
  g1 = parse_vcf(vcf.Reader(filename=args['--vcf']), chrom_list=range(1, 23) + ['X', 'Y'])
