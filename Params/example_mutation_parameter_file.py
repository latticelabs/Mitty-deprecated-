"""Example parameter file for mutate program

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""

chrom = 1  # Chromosome number (VCF file requires this and a simulation might want to use this info)
ref_file = 'Data/adenovirus.smalla'  # Reference sequence in smalla format
out_file = 'Data/variants.vcf'  # VCF file name
models = {
  'snp': {
    'model': 'snp',  # Stock SNP model
    'start': 1000,  # Where on the seq to start generating SNPs
    'stop': 5000,  # Where on the seq to end generating SNPs
    'args': {
      'p': 0.001,  # Stock SNP model has only one parameter
      'poisson_rng_seed': 1,
      'base_sub_rng_seed': 1  # Stock SNP model has two rngs that it needs seeds for
    }
  }
}
