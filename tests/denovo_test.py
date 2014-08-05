from mitty.variation import Variation
from mitty.denovo import *


def arbitrate_variant_collisions_test():
  mask = {c: sparse.lil_matrix((2, 1000), dtype='i1') for c in [1,2]}
  g1 = {1: [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)],
        2: [Variation(7, 10, 'CAA', 'C', HET2)]}  # This should be placed with no problems
  g2 = {1: [Variation(1, 2, 'C', 'CAA', HET2)],  # This will collide
        2: [Variation(7, 8, 'G', 'T', HET1),  # This will pass
            Variation(17, 18, 'G', 'T', HET1)]}  # This will pass
  g1_ = arbitrate_variant_collisions(mask, g1)
  assert g1 == g1_, g1_

  g2_ = arbitrate_variant_collisions(mask, g2)
  assert {1: [], 2: [Variation(7, 8, 'G', 'T', HET1), Variation(17, 18, 'G', 'T', HET1)]}