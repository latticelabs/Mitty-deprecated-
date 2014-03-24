"""This module contains functions that generate variants with reference to a reference genome. The module can be called
as a script as well. This is useful for creating test data for MGR algorithms/data formats. The script will output
fasta file(s) with the mutated chromosome(s) and VFC file(s).

Usage:
mutate [snp] [options] [verbose]

Options:
  --ref=<REF>             Reference sequence [default: porcine_circovirus.fa]
  --psnp=<SNP>            Probability of SNPs [default: 0.01]
  --start=<START>         Where to start on the sequence [default: 0]
  --stop=<STOP>           Where to end on the sequence (-1 means end of the sequence) [default: -1]
  --out=<OUT>             Output file name [default: mutated.fa]
  --seed=<SEED>           Seed for RNG [default: 1]

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
__version__ = '0.1.0'
vcf_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

import seqio, numpy, docopt, logging
logger = logging.getLogger(__name__)

from Bio import SeqIO


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
  1. All mark coordinates (copy end and jump) refer to the reference sequence
  2. At the first command (starts the sequence) the copy end marker is ignored and the jump marker is 0
  3. At the last command (ends the sequence) the jump marker is ignored

  Algorithm

  All the mutations we perform on a sequence can be expressed as combinations of SNP, DEL or INS. Given our mutation
  requirements we generate a 'mutation program' which consists of sequentially ordered markers along the reference
  sequence. The markers are used to indicate which parts of the original sequence are copied to generate the
  new sequence incrementally.


  0                50 51        100         200           300              500
  |-----------------|-|----------|-----------|-------------|----------------|
  |                 | |          |           |             |                |
 start              SNP              DEL                  INS              end


  We start by copying over 0->49
  We insert SNP at 50
  We skip over to 51 and continue copying 51->99
  We skip over to 200 and continue copying 200->300
  At 300 we insert a sequence and then continue copying 300->499 (end)

  We carry with us a pair of start and stop pointers which tell us what to copy next, while each pointer carries
  information on what to do inbetween the copying.

  SNP - we have a single base that we insert before continuing copying
  DEL - we just skip over
  INS - we carry a sequence that we insert here before continuing copying

  """
  mut_seq = bytearray()
  for n in range(len(mutation_program)-1):
    mut_seq += ref_seq[mutation_program[n][1]:mutation_program[n+1][0]]
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
  lam = 1./p  # The expected interval between SNPs is the inverse of the probability that a SNP is going to occur
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


def create_mutation_program(ref_seq, snp_commands):
  """See notes from 'polymerize'."""
  mut_prg = [(None, 0, None)]  # Start command
  if snp_commands is not None:
    for snp_c in snp_commands:
      mut_prg.append((snp_c[0], snp_c[0]+1, snp_c[1]))
  mut_prg.append((len(ref_seq), None, None))  # End command
  return mut_prg


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
    #CHROM POS     ID        REF    ALT     QUAL FILTER INFO""".format(sim_date, __version__, reference_filename)
  )

def write_vcf_mutations(file_handle, mutations):
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
  for n in range(len(mutations['CHROM'])):
    file_handle.write("\t".join([mutations[k] for k in vcf_columns]))

def main(args):
  reference = SeqIO.parse(args['-r'],'fasta').next() #We always have a reference
  if args['snp']:
    parameters = {
      'p': float(args['-p'])/100,  # Convert from percentage to fraction
      'sm': None  # TODO
    }
    block_size = int(args['--block_size'])
    snp(reference, parameters, block_size)

if __name__ == "__main__":
  args = docopt.docopt(__doc__, version=__version__)
  level = logging.DEBUG if args['verbose'] else logging.WARNING
  logging.basicConfig(level=level)

  with open(args['--ref'],'r') as f:
    header, ref_seq = seqio.fast_read_fasta(f)

  # TODO Push all parameters to config file
  if args['snp']:
    snp_commands = create_snps(ref_seq, float(args['--psnp']), block_size=100000, seed=int(args['--seed']))
  else:
    snp_commands = None

  mutation_program = create_mutation_program(ref_seq, snp_commands)
  mutated_seq = polymerize(ref_seq, mutation_program)
  mutated_header = header + ' Mutated by mutate {:s}'.format(__version__)

  with open(args['--out'], 'w') as f:
    seqio.fast_write_fasta(f, mutated_seq, mutated_header, width=70)