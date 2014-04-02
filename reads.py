"""This module contains functions that generate simulated reads. The module can be called
as a script as well. This is useful for creating test data for MGR algorithms/data formats

Usage:
reads [options] [verbose]
reads formats

Options:
  --ref=REF                  Reference sequence [default: Data/porcine_circovirus.fa]
  --vcf=VCF                  VCF file (If none, null model is generated)
  --bam=BAM                  Output BAM file [default: Data/test.bam]
  --seed=SEED                Seed for random number generator [default: 0]
  --read_profile=RP          A profile file that describes characteristics for our reads
  --block_len=BL             Block length (adjust according to compute resources) [default: 1000000]
  --comments=COM             User comments to be saved in side car file
  verbose                    Dump detailed logger messages

Notes:
1. The quality scores are in Phred scale (as specified in the SAM spec)
2. Any annotations the user wishes to make (plus the command line arguments and all other parameters used to run the
   sim) are stored in a sidecar file with the same name as the bam file with .info added to the end

"""
#__version__ = '0.1.0'  # Fixed read lengths. Uniform coverage
#__version__ = '0.2.0'  # Paired end. Fixed read lengths. Uniform coverage
__version__ = '0.2.1'  # Paired end. Fixed read lengths. Uniform coverage. Handles VCF files, does not regenerate whole
# sequences. Ignores multiple variants

import imp
import docopt
import numpy
from Bio import SeqIO  # Needed for loading reference sequence
import vcf  # Needed for handling VCF files
import pysam  # Needed to write BAM files
import logging

logger = logging.getLogger(__name__)


def main(args):
  """
  1. Prepare the files for access
  2. Load parameters
  3. For each block
     a) Load the sequence
     b) Load the relevant variants
     c) Generate the reads and write to BAM file
  4. Write the sidecar file
  5. Cleanup and exist
  """
  rng = numpy.random.RandomState(int(args['--seed']))

  # Load the read profile as a module
  params = imp.load_source('params', args['--read_profile'], open(args['--read_profile'], 'r'))

  ref = SeqIO.read(args['--ref'], 'fasta')
  ref_description = ref.description
  ref_seq = ref.seq.tostring()
  # Load the sequence
  # TODO: Load in blocks? Needed?
  # TODO: Error checking?

  vcf_reader = vcf.Reader(filename=args['--vcf']) if args['--vcf'] is not None else None
  chrom = 1

  bam_hdr = {'HD': {'VN': '1.4'},
             'SQ': [{'LN': 300000000, 'SN': ref.description, 'SP': 'simulated human'}]}
  bam_file = pysam.Samfile(args['--bam'], 'wb', header=bam_hdr)  # Write binary BAM with header
  ref_len = len(ref_seq)
  blk_start = 0
  blk_len = int(args['--block_len'])
  while blk_start < ref_len:
    blk_stop = min(ref_len, blk_start + blk_len)
    variants = {rec.POS: rec for rec in
                vcf_reader.fetch(chrom, blk_start, blk_stop)} if vcf_reader is not None else {}
    # Fetch all variants in this block and structure as a dictionary keyed by coordinate
    seqs = polymerize(ref_seq[blk_start:blk_stop], blk_start, variants)  # Generate the sequences for this block
    reads = generate_reads(
      seqs=seqs,
      read_len=params.read_len,
      template_len=params.template_len,
      coverage=params.coverage,
      paired=params.paired,
      rng=rng)  # Generate the reads from this block
    # List of tuples. Paired reads will have two reads per tuple, otherwise only one
    save_reads_to_bam(bam_file, reads, params.template_len)
    blk_start += blk_len
  bam_file.close()

  with open(args['--bam'] + '.info', 'w') as file_handle:
    write_sidecar(args, args['--read_profile'], file_handle)


def polymerize(ref_seq_block, start_coord, variant_dict):
  """Given a part of the reference sequence and variants, generate as many variant sequences as called for
  """
  # Place holder: copy the sequence element by element testing if a variant exists there and then implementing it
  # Only does SNPs
  mutated_seq_block = bytearray()
  # TODO: Make this calculation more efficient
  for n in range(len(ref_seq_block)):
    if start_coord + n in variant_dict:
      mutated_seq_block += variant_dict[start_coord + n].ALT[0].sequence
      if len(variant_dict[start_coord + n].ALT) > 1:
        logger.debug('reads.py does not currently handle multiple variants. Only the first variant is processed.')
    else:
      mutated_seq_block += ref_seq_block[n]

  return [mutated_seq_block.__str__()]


def generate_reads(seqs=[],
                   read_len=None,
                   template_len=None,
                   coverage=5,
                   paired=False,
                   rng=None):
  """Given a list of sequences generate reads with the given characteristics

  Inputs:
    seqs             - list of string(like)s containing the variant DNA sequences
    read_len         - Fixed read length
    template_len     - Template length. Only needed if paired is True
    coverage         - the coverage that we want
    paired           - paired reads or not
    rng              - numpy.random.RandomState(seed)

  Outputs
                                 _________ ( seq_str, quality_str, coordinate)
    reads     -  [              /
                  [( ... ), ( ...)],
                  [( ... ), ( ...)], -> inner list = 2 elements if paired reads, 1 otherwise
                       .
                       .
                       .
                 ] -> outer list = number of reads

  Quality: Sanger scale 33-126
  """
  logger.debug('Starting to generate reads')
  # Placeholder, generate reads from first sequence only
  seq = seqs[0]  #Sequence being considered
  num_reads = int(coverage * (float(len(seq)) / float(read_len)))
  rl = read_len
  tl = template_len
  if paired:
    rd_st = rng.randint(0, len(seq) - template_len, num_reads)  # Reads are uniformly distributed
    reads = [[(seq[rd_st[n]:rd_st[n] + rl], '~' * rl, rd_st[n]),
              (seq[rd_st[n] + tl - rl:rd_st[n] + tl], '~' * rl, rd_st[n] + tl - rl)] for n in range(num_reads)]
  else:
    rd_st = rng.randint(0, len(seq) - read_len, num_reads)
    reads = [[(seq[rd_st[n]:rd_st[n] + rl], '~' * rl, rd_st[n])] for n in range(num_reads)]
  logger.debug('Finished generating reads')
  return reads


def corrupt_reads(reads):
  """Placeholder, does not damage reads yet. Should corrupt reads and update the quality string. Should get the
  corruption program externally"""
  return reads


def save_reads_to_bam(bam_file, reads, template_len):
  """
  Inputs:
                                 _________ ( seq_str, quality_str, coordinate)
    reads     -  [              /
                  [( ... ), ( ...)],
                  [( ... ), ( ...)], -> inner list = 2 elements if paired reads, 1 otherwise
                       .
                       .
                       .
                 ] -> outer list = number of reads
  """
  if len(reads) == 0: return
  paired = True if len(reads[0]) == 2 else False
  for n in range(len(reads)):
    ar = pysam.AlignedRead()
    ar.qname = "r{:d}".format(n)  # Figure out how to store coordinate
    ar.seq = reads[n][0][0]
    ar.qual = reads[n][0][1]
    if paired:
      ar.flag = 0x41  # end1 0x01 flag has to be set to indicate multiple segments
      ar.tlen = template_len
    bam_file.write(ar)

    if paired:
      ar.seq = reads[n][1][0]
      ar.qual = reads[n][1][1]
      ar.flag = 0x81  # end2 0x01 flag has to be set to indicate multiple segments
      ar.tlen = template_len
      bam_file.write(ar)


def write_sidecar(args, params_file, file_handle):
  """
  Write parameters into a sidecar file
  Inputs:
    args        - the program arguments as parsed by docopts
    file_handle - handle of an opened text file. The output will be appended to this file.
  """
  file_handle.write('SBG Read simulator v{:s}\n'.format(__version__))
  file_handle.write('Commandline:\n')
  for k,v in args.iteritems():
    file_handle.write('{:s}: {:s}\n'.format(k, str(v)))
  file_handle.write('Parameters:\n')
  with open(params_file, 'r') as f:
    for line in f.readlines():
      file_handle.write(line)


if __name__ == "__main__":
  arguments = docopt.docopt(__doc__, version=__version__)
  logging.basicConfig(level=logging.DEBUG)
  main(arguments)