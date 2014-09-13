from mitty.denovo import *
from mitty.tests import *


def merge_test1():
  """Merge variants, non overlapping, existing first (ED)."""
  c1 = [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)]
  dnv = [Variation(10, 13, 'CAA', 'C', HOMOZYGOUS)]
  c2 = merge_variants_with_chromosome(c1, dnv)
  assert c2 == [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS), Variation(10, 13, 'CAA', 'C', HOMOZYGOUS)]


def merge_test2():
  """Merge variants, non overlapping, denovo first (DE)."""
  dnv = [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)]
  c1 = [Variation(10, 13, 'CAA', 'C', HOMOZYGOUS)]
  c2 = merge_variants_with_chromosome(c1, dnv)
  assert c2 == [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS), Variation(10, 13, 'CAA', 'C', HOMOZYGOUS)]


def merge_test3():
  """Merge variants, non overlapping (EDED)"""
  c1 = [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS), Variation(13, 16, 'CTT', 'C', HOMOZYGOUS)]
  dnv = [Variation(8, 11, 'CCC', 'C', HOMOZYGOUS), Variation(20, 23, 'CAA', 'C', HOMOZYGOUS)]
  c2 = merge_variants_with_chromosome(c1, dnv)
  assert c2 == [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS),
                Variation(8, 11, 'CCC', 'C', HOMOZYGOUS),
                Variation(13, 16, 'CTT', 'C', HOMOZYGOUS),
                Variation(20, 23, 'CAA', 'C', HOMOZYGOUS),]


def merge_test4a():
  """Merge variants, full overlapping (E-D-)"""
  c1 = [Variation(2, 5, 'CAA', 'C', HOMOZYGOUS)]
  dnv = [Variation(2, 5, 'CAA', 'T', HOMOZYGOUS)]
  c2 = merge_variants_with_chromosome(c1, dnv)
  assert c2 == [Variation(2, 5, 'CAA', 'C', HOMOZYGOUS)]


def merge_test4b():
  """Merge variants, full overlapping, SNP (D-E-)."""
  c1 = [Variation(2, 3, 'C', 'G', HET1)]
  dnv = [Variation(2, 3, 'C', 'T', HET1)]
  c2 = merge_variants_with_chromosome(c1, dnv)
  assert c2 == [Variation(2, 3, 'C', 'G', HET1)]


# Important test - killed a nasty logic bug
def merge_test4c():
  """Merge variants, full overlapping, with non-colliding preceder"""
  c1 = [Variation(1, 2, 'G', 'C', HET2),
        Variation(2, 3, 'A', 'C', HET1)]
  dnv = [Variation(2, 3, 'A', 'T', HET1)]
  correct_final_c = [Variation(1, 2, 'G', 'C', HET2),
                     Variation(2, 3, 'A', 'C', HET1)]
  c2 = merge_variants_with_chromosome(c1, dnv)
  assert c2 == correct_final_c, c2


def merge_test4():
  """Merge variants, overlapping (E-D). D will collide and will be rejected"""
  c1 = [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)]
  dnv = [Variation(2, 5, 'CCC', 'C', HOMOZYGOUS)]
  c2 = merge_variants_with_chromosome(c1, dnv)
  assert c2 == [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)]


def merge_test5():
  """Merge variants, overlapping (D-E). D will collide and will be rejected"""
  dnv = [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)]
  c1 = [Variation(2, 5, 'CCC', 'C', HOMOZYGOUS)]
  c2 = merge_variants_with_chromosome(c1, dnv)
  assert c2 == [Variation(2, 5, 'CCC', 'C', HOMOZYGOUS)]


def merge_test6():
  """Merge variants, overlapping heterozygous. Overlapping but with mixed zygosity"""
  c1 = [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS),
        Variation(13, 16, 'CTT', 'C', HET1),
        Variation(20, 23, 'CAA', 'C', HET1),
        Variation(26, 29, 'CGG', 'C', HOMOZYGOUS)]
  dnv = [Variation(5, 8, 'CCC', 'C', HOMOZYGOUS),  # Will collide out
         Variation(13, 16, 'CAA', 'C', HET2),  # Will pass
         Variation(20, 23, 'CAA', 'C', HOMOZYGOUS),  # Will collide out
         Variation(26, 29, 'CCC', 'C', HET2)]  # Will collide out
  c2 = merge_variants_with_chromosome(c1, dnv)
  assert c2 == [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS),
                Variation(13, 16, 'CTT', 'C', HET1),
                Variation(13, 16, 'CAA', 'C', HET2),
                Variation(20, 23, 'CAA', 'C', HET1),
                Variation(26, 29, 'CGG', 'C', HOMOZYGOUS)], c2


def add_variant_model_to_genome_test1():
  g1 = {1: [Variation(12, 13, 'G', 'C', HET2),
            Variation(13, 14, 'A', 'C', HET1),
            Variation(19, 20, 'G', 'A', HET1)]}  # This should be placed with no problems

  def variant_generator():
    g2 = [{1: [Variation(13, 14, 'A', 'T', HET1)]}]  # This will pass
    for g in g2:
      yield g

  correct_final_g = {1: [Variation(12, 13, 'G', 'C', HET2),
                         Variation(13, 14, 'A', 'C', HET1),
                         Variation(19, 20, 'G', 'A', HET1)]}

  vg = variant_generator()
  g2 = merge_variants_with_genome(g1, vg)

  assert correct_final_g == g2, g2


def add_variant_model_to_genome_test1a():
  g1 = {1: [Variation(13, 14, 'A', 'C', HET1)]}  # This should be placed with no problems

  def variant_generator():
    g2 = [{1: [Variation(13, 14, 'A', 'T', HET1)]}]  # This will pass
    for g in g2:
      yield g

  correct_final_g = {1: [Variation(13, 14, 'A', 'C', HET1)]}

  vg = variant_generator()
  g2 = merge_variants_with_genome(g1, vg)

  assert correct_final_g == g2, g2


def add_variant_model_to_genome_test2():
  g1 = {1: [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)],
        2: [Variation(7, 10, 'CAA', 'C', HET2)]}  # This should be placed with no problems

  def variant_generator():
    g2 = [{1: [Variation(1, 2, 'C', 'CAA', HET2)]},  # This will collide
          {2: [Variation(7, 8, 'G', 'T', HET1)]},  # This will pass
          {2: [Variation(17, 18, 'G', 'T', HET1)]}]  # This will pass
    for g in g2:
      yield g

  correct_final_g = {1: [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)],
                     2: [Variation(7, 10, 'CAA', 'C', HET2),
                         Variation(7, 8, 'G', 'T', HET1),
                         Variation(17, 18, 'G', 'T', HET1)]}

  vg = variant_generator()
  g2 = merge_variants_with_genome(g1, vg)

  assert correct_final_g == g2, g2
  assert g1 == {1: [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)], 2: [Variation(7, 10, 'CAA', 'C', HET2)]}  # Must not change original


# This test uses the stock SNP which must exist for this test to pass
def load_variant_models_test():
  """Loading SNP variant as a test."""
  param_json = [
        {
          "snp": {
             "phet": 0.5,
             "p": 0.01,
             "master_seed": 1
          }
        },
        {
          "snp": {
             "phet": 0.0,
             "p": 0.01,
             "master_seed": 2
          }
        }
      ]
  mdl = load_variant_models(param_json)
  assert hasattr(mdl[0]["model"], 'variant_generator')
  assert 'master_seed' in mdl[1]["params"]


# This test uses the stock SNP which must be functional for this test to pass as well as some files created during
# test setup
def main_test():
  """Integrating denovo and stock SNP"""
  import json, vcf, tempfile
  param_json = {
    "denovo_variant_models": [
      {
        "snp": {
           "chromosome": [1, 2],
           "phet": 0.5,
           "p": 0.0001,
           "master_seed": 1
        }
      },
      {
        "snp": {
           "chromosome": [2],
           "phet": 0.0,
           "p": 0.01,
           "master_seed": 2
        }
      }
    ]
  }
  _, vcf_file_fname = tempfile.mkstemp(suffix='.vcf.gz')
  ref = FastaGenome(seq_dir=example_fasta_genome)
  g1 = main(ref, vcf_file_name=vcf_file_fname, parameters=param_json, master_seed=1)
  assert os.path.exists(vcf_file_fname)

  vcf_rdr = vcf.Reader(filename=vcf_file_fname).fetch(chrom=1, start=0)
  v1 = vcf_rdr.next()
  assert v1.POS, g1[1][0].POS