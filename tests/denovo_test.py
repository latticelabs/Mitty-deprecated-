from mitty.variation import Variation
from mitty.denovo import *
from . import *


def arbitrate_variant_collisions_test():
  mask = init_mask({c: 1000 for c in [1, 2]})
  g1 = {1: [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)],
        2: [Variation(7, 10, 'CAA', 'C', HET2)]}  # This should be placed with no problems
  g2 = {1: [Variation(1, 2, 'C', 'CAA', HET2)],  # This will collide
        2: [Variation(7, 8, 'G', 'T', HET1),  # This will pass
            Variation(17, 18, 'G', 'T', HET1)]}  # This will pass
  g1_ = arbitrate_variant_collisions(g1, mask)
  assert g1 == g1_, g1_

  g2_ = arbitrate_variant_collisions(g2, mask)
  assert {1: [], 2: [Variation(7, 8, 'G', 'T', HET1), Variation(17, 18, 'G', 'T', HET1)]} == g2_, g2_


def add_variants_to_genome_test():
  mask = init_mask({c: 1000 for c in [1, 2]})
  g1 = {1: [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)],
        2: [Variation(7, 10, 'CAA', 'C', HET2)]}  # This should be placed with no problems
  fill_mask(mask, g1)

  def variant_generator():
    g2 = [{1: [Variation(1, 2, 'C', 'CAA', HET2)]},  # This will collide
          {2: [Variation(7, 8, 'G', 'T', HET1)]},  # This will pass
          {2: [Variation(17, 18, 'G', 'T', HET1)]}]  # This will pass
    for g in g2:
      yield g

  correct_final_g = {1: [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)],
                     2: [Variation(7, 10, 'CAA', 'C', HET2), Variation(7, 8, 'G', 'T', HET1),
                         Variation(17, 18, 'G', 'T', HET1)]}

  vg = variant_generator()
  add_variants_to_genome(g1, mask, vg)

  assert correct_final_g == g1, g1


# This test uses the stock SNP which must exist for this test to pass
def load_variant_models_test():
  """Loading SNP variant as a test."""
  param_json = {
    "variant_models": [
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
  }
  mdl = load_variant_models(param_json)
  assert hasattr(mdl[0]["model"], 'variant_generator')
  assert 'master_seed' in mdl[1]["params"]


# This test uses the stock SNP which must be functional for this test to pass as well as some files created during
# test setup
def main_test():
  import json, vcf
  param_json = {
    "variant_models": [
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
  _, param_fname = tempfile.mkstemp(suffix='.json')
  with open(param_fname, 'w') as fp:
    json.dump(param_json, fp, indent=2)
  _, vcf_file_fname = tempfile.mkstemp(suffix='.vcf.gz')
  g1 = main(wg_name, vcf_file_name=vcf_file_fname, param_file_name=param_fname, master_seed=1)
  assert os.path.exists(vcf_file_fname)

  vcf_rdr = vcf.Reader(filename=vcf_file_fname).fetch(chrom=1, start=0)
  v1 = vcf_rdr.next()
  assert v1.POS, g1[1][0].POS
