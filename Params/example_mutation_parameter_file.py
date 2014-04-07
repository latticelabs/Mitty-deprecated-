"""Example parameter file for mutate program

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""

chrom = '1'  # Chromosome number (VCF file requires this and a simulation might want to use this info)
ref = 'Data/porcine_circovirus.smalla'  # Reference sequence in smalla format
vcf = 'Data/variants.vcf'  # VCF file out (Set to None to write to stdout)
models = {
  'snp': {
    'model': 'snp',  # Stock SNP model
    'start': 0,  # Where on the seq to start generating SNPs
    'stop': -1,  # Where on the seq to end generating SNPs
    'args': {
      'p': 0.01,  # Stock SNP model has only one parameter
      'poisson_rng_seed': 1,
      'base_sub_rng_seed': 1  # Stock SNP model has two rngs that it needs seeds for
    }
  }
}
