"""This module contains functions that generate variants with reference to a reference genome. The module can be called
as a script as well. This is useful for creating test data for MGR algorithms/data formats

Usage:
mutate snp [-r <REF>] [-p <SNP>] [--block_size <BS>]

Options:
  -r=<REF>                Reference chromosome [default: ../../Data/GRCh38/chr24.fa]
  -p=<SNP>                What percentage of bases (on average) are mutated [default: 0.5]
  --block_size==<BS>      Computation block size [default: 1000000]

The mutations are carried internally ("mutator format") as a dictionary version of the VCF format. The dictionary
contains keys corresponding to the VCF headers and each value is a list of values corresponding to each mutation we
simulate.
"""

import numpy, docopt
from Bio import SeqIO

def snp(reference, parameters, block_size=1000000):
  """
  Inputs:
    reference  - reference sequence string
    parameters - dictionary
     'p': the rate of mutations.
         0 = no mutation 1 = every base is mutated
     'sm' : substitution matrix.
         4x4 list of lists indicating relative probabilities of one base mutating into another
         Rows are the original bases, columns are the targets. Each row (disregarding the diagonal) must add up to 1.0
         Order of indexes is ATCG.
         If None then all mutations are equi-probable.
    block_size  - how many bases do we process at a time. Related to our system resources.

  Outputs:
    Variant calls in the mutator format

  Algorithm:
    The base rate of mutations is p. We could go through each base position, toss a biased coin to decide if that
    position should be mutated. It is more efficient to note that this is a poisson point process and to simply draw
    numbers from a canned poisson generator in blocks to determine which positions are to be mutated.
  """

def main(args):
  reference = SeqIO.parse(args['-r'],'fasta').next() #We always have a reference
  if args['snp']:
    parameters = {
      'p': float(args['-p']),
      'sm': None  # TODO
    }
    block_size = int(args['--block_size'])
    snp(reference, parameters, block_size)

if __name__ == "__main__":
  arguments = docopt.docopt(__doc__, version='v1')
  main(arguments)