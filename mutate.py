"""This module contains functions that generate variants with reference to a reference genome. The module can be called
as a script as well. This is useful for creating test data for MGR algorithms/data formats. The script will output
fasta file(s) with the mutated chromosome(s) and VCF file(s).

Usage:
mutate [snp] [options] [verbose]

Options:
  snp                     Generate SNPs
  --ref=REF               Reference sequence [default: porcine_circovirus.fa]
  --out=OUT               Output file name [default: mutated]
  --start=START           Where to start on the sequence [default: 0]
  --stop=STOP             Where to end on the sequence (-1 means end of the sequence) [default: -1]
  --seed=SEED             Seed for RNG [default: 1]
  --block_size=BS         Block size for operations. Adjust to match memory/resources of platform [default: 100000]
  --paramfile=PFILE       Name for parameter file [default: Params/example_mutation_parameter_file.py]
  verbose                 Dump detailed logger messages

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
__version__ = '0.1.0'

import seqio
import numpy
import docopt
import datetime
import logging

logger = logging.getLogger(__name__)


def polymerize(ref_seq, mutation_program):
  """This function takes in a sequence and a mutation program and generates a mutated sequence
  Inputs:
    ref_seq           - reference sequence
    mutation_program  - (see below)

  Output:
    mut_seq           - mutated sequence

  The mutation program is a list of tuples representing commands:
    The first element is the copy end marker - where should we stop copying at
    The second element is the jump marker    - where should we carry on copying from
    The third element is the sequence that should be inserted before we continue copying. This is None for 'dels'

  Notes:
  0. The program must be in ascending order
  1. All mark coordinates (copy end and jump) refer to the reference sequence
  2. At the first command (starts the sequence) the copy end marker is ignored and the jump marker is 0
  3. At the last command (ends the sequence) the jump marker is ignored

  Algorithm

  All the mutations we perform on a sequence can be expressed as combinations of DEL and INS. Given our mutation
  requirements we generate a 'mutation program' which consists of sequentially ordered markers along the reference
  sequence. The markers are used to indicate which parts of the original sequence are copied to generate the
  new sequence incrementally and what sequences are inserted where.

  Say we have a 500bp sequence and we have the following mutations

  SNP @10
  Deletion 25-75
  Insertion @100
  Tandem repeat of 125-150
  Interspersed repeat of 175-200 at 225
  Translocation from 250-275 to 300

  These mutations can be represented as insertions and deletions as follows

                   SNP           deletion       Insertion      Tandem repeat     Interspersed      Translocation ...

  0              9    11        24      75      100   101      150   151           225   226      249  275      300
  |--------------|     |--------|       |--------|     |--------|     |-------------|     |--------|    |--------|
  |    copy      |.....|  copy  |       |  copy  |.....|  copy  |.....|    copy     |.....|  copy  |    |  copy  |....
                   INS                             INS            INS                 INS                          INS

  For every mutation event we carry a stop copy and resume copy pointers and a sequence (if an insertion is called for)
  At the beginning and end we add two dummy events to ensure the beginning and end of the sequence are properly copied.

  We copy from the resume copy of event n to the stop copy of event n+1
  If there is a insertion sequence we insert that before repeating from n+1 and so on
  """
  mut_seq = bytearray()
  for n in range(1, len(mutation_program)):
    mut_seq += ref_seq[mutation_program[n - 1][1]:mutation_program[n][0]]
    if mutation_program[n][2] is not None: mut_seq += mutation_program[n][2]
  return mut_seq


def create_snps(ref_seq, p, block_size=100000, seed=1):
  """Given a reference sequence and a snp probability generate snp commands.
  Inputs:
    ref_seq      - reference sequence
    p            - probability of any one base becoming mutated
    block_size   - we generate these many SNPs at a time until we are done with the sequence
    seed         - seed for RNG

  Output:
    snp_commands - list of tuples:
                   First element is location of SNP
                   Second element is substituted base

  Algorithm:
    The base rate of mutations is p. We could go through each base position, toss a biased coin to decide if that
    position should be mutated. It is more efficient to note that this is a poisson point process and to simply draw
    numbers from a canned poisson generator in blocks to determine which positions are to be mutated.

  Future expansion:
  1. Non uniform SNPs
  2. Non uniform substitution matrix
  """
  base_sub_mat = {
    'G': 'ATC',
    'A': 'TCG',
    'T': 'CGA',
    'C': 'GAT'
  }
  rng = numpy.random.RandomState(seed)  # Initialize the numpy RNG
  lam = 1. / p  # The expected interval between SNPs is the inverse of the probability that a SNP is going to occur
  snp_loc = 0
  snp_commands = []
  while snp_loc < len(ref_seq):
    snp_intervals = rng.poisson(lam=lam, size=block_size)
    base_subs = rng.randint(3, size=block_size)
    for snpi, bsub in zip(snp_intervals, base_subs):
      snp_loc += snpi
      if snp_loc >= len(ref_seq): break
      snp_commands.append((snp_loc, base_sub_mat[chr(ref_seq[snp_loc])][bsub]))
  logger.debug('{:d} SNPs'.format(len(snp_commands)))
  return snp_commands


def resolve_conflicts_and_create_mutation_program(ref_seq, snp_commands):
  """Transcribe the SNP and SV commands into a mutation program and a list of vcf strings.
  Inputs:
    ref_seq      - reference sequence
    snp_commands - list of tuples:
                   First element is location of SNP
                   Second element is substituted base
  Output:
    mutation_program  - The mutation program is a list of tuples representing commands:
          The first element is the copy end marker - where should we stop copying at
          The second element is the jump marker    - where should we carry on copying from
          The third element is the sequence that should be inserted before we continue copying. This is None for 'dels'
    vcf_lines     - Lines that correspond to information needed for printing in the VCF format


  See notes from 'polymerize'.
  Right now this function simply transcribes SNP commands into mutation program. When we have SVs too we will do
  something slightly more involved where we a) resolve conflicts between mutations and b) make sure the program is
  created in ascending coordinate order.
  """
  mut_prg = [(None, 0, None)]  # Start command
  vcf_lines = []
  if snp_commands is not None:
    for snp_c in snp_commands:
      mut_prg.append((snp_c[0], snp_c[0] + 1, snp_c[1]))
      vcf_lines.append(['1',  # CHROM We default the chrom no to 1
                        str(snp_c[0] + 1),  # POS The numbering for bases starts at 1
                        '.',  # No id - this is fake
                        chr(ref_seq[snp_c[0]]),  #REF
                        snp_c[1],  # ALT
                        '96',  # Arbitrary, high, Phred score
                        'PASS',  # Passed the filters - fake mutation
                        '.'])  # Read depth unknown
  mut_prg.append((len(ref_seq), None, None))  # End command
  return mut_prg, vcf_lines


def write_vcf_header(file_handle, sim_date, reference_filename):
  """Given a file handle, write out a suitable header to start the VCF file
  Inputs:
    file_handle       - an open file handle

  Notes:
  1. 'filedate' is the date of the simulation
  2. 'source' contains the version of the mutate program used to generate the data
  3. 'reference' contains the name of the file entered as the reference file

  """
  file_handle.write(
    """##fileformat=VCFv4.1
##fileDate={:s}
##source=mutate {:s}
##reference={:s}
#CHROM POS     ID        REF    ALT     QUAL FILTER INFO\n""".format(sim_date, __version__, reference_filename)
  )


def write_vcf_mutations(file_handle, vcf_lines):
  """Given a mutator format dictionary write the mutations in VCF format into the file
  Inputs:
    file_handle   - handle of an opened text file. The output will be appended to this file.
    mutations     - dictionary of lists in mutator format
  Notes:
  1. The mutator format is a dictionary of lists which makes it easy to drop the information into the VCF file.
     The dictionary keys correspond to the VCF columns. Each value is a list with however many elements we want the
     VCF file to have
  2. The ID is a madeup id that is guaranteed not to clash with any dbSNP id because it is a string unique to the
     mutate program
  """
  for line in vcf_lines:
    file_handle.write("\t".join(line) + '\n')


if __name__ == "__main__":
  args = docopt.docopt(__doc__, version=__version__)
  level = logging.DEBUG if args['verbose'] else logging.WARNING
  logging.basicConfig(level=level)
  logging.debug(args)

  pars = {}
  execfile(args['--paramfile'], {}, pars)

  with open(args['--ref'], 'r') as f:
    header, ref_seq = seqio.fast_read_fasta(f)

  if pars.has_key('snp'):
    snp_commands = create_snps(ref_seq, pars['snp']['p'], block_size=int(args['--block_size']),
                               seed=int(args['--seed']))
  else:
    snp_commands = None

  mutation_prog, vcf_lines = resolve_conflicts_and_create_mutation_program(ref_seq, snp_commands)
  mutated_seq = polymerize(ref_seq, mutation_prog)
  mutated_header = 'Mutated by mutate {:s} '.format(__version__) + header

  with open(args['--out'] + '.fa', 'w') as f:
    seqio.fast_write_fasta(f, mutated_seq, mutated_header, width=70)

  with open(args['--out'] + '_params.txt', 'w') as f:
    f.write('Commandline parameters: \n' + args.__str__() + '\n\nParameter file: \n' + pars.__str__())

  with open(args['--out'] + '_variants.vcf', 'w') as f:
    write_vcf_header(f, datetime.datetime.now().isoformat(), args['--ref'])
    write_vcf_mutations(f, vcf_lines)