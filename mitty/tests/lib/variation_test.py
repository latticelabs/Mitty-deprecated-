import mitty.lib.variation as vr
from mitty.lib.variation import HOMOZYGOUS, HET_01, HET_10
from nose.tools import assert_sequence_equal

nv = vr.new_variation


def copy_test1():
  """Copy list of variants with no modification"""
  c1 = [v0, v1] = [nv(1, 4, 'CAA', 'C', HOMOZYGOUS), nv(13, 16, 'CTT', 'C', HOMOZYGOUS)]
  c2 = vr.copy_variant_sequence(c1)
  assert_sequence_equal(c1, list(c2))


def copy_test2():
  """Copy list of variants selectively"""
  c1 = [v0, v1, v2] = [nv(1, 4, 'CAA', 'C', HOMOZYGOUS),
                       nv(7, 8, 'C', 'T', HOMOZYGOUS),
                       nv(13, 16, 'CTT', 'C', HOMOZYGOUS)]
  c2 = vr.copy_variant_sequence(c1, idx=[0, 2])
  assert_sequence_equal([v0, v2], list(c2))


def copy_test3():
  """Copy list of variants selectively with zygosity"""
  c1 = [v0, v1, v2] = [nv(1, 4, 'CAA', 'C', HOMOZYGOUS),
                       nv(7, 8, 'C', 'T', HOMOZYGOUS),
                       nv(13, 16, 'CTT', 'C', HOMOZYGOUS)]
  c2 = vr.copy_variant_sequence(c1, idx=[0, 2], het=[HET_01, HET_10])
  v0d = vr.Variation(v0.vd, het=HET_01)
  v2d = vr.Variation(v2.vd, het=HET_10)
  assert_sequence_equal([v0d, v2d], list(c2))


def merge_test1():
  """Merge variants, non overlapping, existing first (ED)."""
  c1 = [v10] = [nv(1, 4, 'CAA', 'C', HOMOZYGOUS)]
  c2 = [v20] = [nv(10, 13, 'CAA', 'C', HOMOZYGOUS)]
  c3 = vr.merge_variants(c1, c2)
  assert_sequence_equal(c3, [v10, v20])


def merge_test2():
  """Merge variants, non overlapping, denovo first (DE)."""
  c1 = [v10] = [nv(10, 13, 'CAA', 'C', HOMOZYGOUS)]
  c2 = [v20] = [nv(1, 4, 'CAA', 'C', HOMOZYGOUS)]
  c3 = vr.merge_variants(c1, c2)
  assert_sequence_equal(c3, [v20, v10])


def merge_test3():
  """Merge variants, non overlapping (EDED)"""
  c1 = [v10, v11] = [nv(1, 4, 'CAA', 'C', HOMOZYGOUS), nv(13, 16, 'CTT', 'C', HOMOZYGOUS)]
  c2 = [v20, v21] = [nv(8, 11, 'CCC', 'C', HOMOZYGOUS), nv(20, 23, 'CAA', 'C', HOMOZYGOUS)]
  c3 = vr.merge_variants(c1, c2)
  assert_sequence_equal(c3, [v10, v20, v11, v21])


def merge_test4a():
  """Merge variants, full overlapping (E-D-)"""
  c1 = [v10] = [nv(2, 5, 'CAA', 'C', HOMOZYGOUS)]
  c2 = [v20] = [nv(2, 5, 'CAA', 'T', HOMOZYGOUS)]
  c3 = vr.merge_variants(c1, c2)
  assert_sequence_equal(c3, [v10])


def merge_test4b():
  """Merge variants, full overlapping, SNP (D-E-)."""
  c1 = [v10] = [nv(2, 3, 'C', 'G', HET_10)]
  c2 = [v20] = [nv(2, 3, 'C', 'T', HET_10)]
  c3 = vr.merge_variants(c1, c2)
  assert_sequence_equal(c3, [v10])


# Important test - killed a nasty logic bug
def merge_test4c():
  """Merge variants, full overlapping, with non-colliding preceder"""
  c1 = [v10, v11] = [nv(1, 2, 'G', 'C', HET_01),
                     nv(2, 3, 'A', 'C', HET_10)]
  c2 = [v20] = [nv(2, 3, 'A', 'T', HET_10)]
  c3 = vr.merge_variants(c1, c2)
  assert_sequence_equal(c3, [v10, v11])


def merge_test4():
  """Merge variants, overlapping (E-D). D will collide and will be rejected"""
  c1 = [v10] = [nv(1, 4, 'CAA', 'C', HOMOZYGOUS)]
  c2 = [v20] = [nv(2, 5, 'CCC', 'C', HOMOZYGOUS)]
  c3 = vr.merge_variants(c1, c2)
  assert_sequence_equal(c3, [v10])


def merge_test5():
  """Merge variants, overlapping (D-E). D will collide and will be rejected"""
  c2 = [v20] = [nv(1, 4, 'CAA', 'C', HOMOZYGOUS)]
  c1 = [v10] = [nv(2, 5, 'CCC', 'C', HOMOZYGOUS)]
  c3 = vr.merge_variants(c1, c2)
  assert_sequence_equal(c3, [v10])


def merge_test6():
  """Merge variants, overlapping heterozygous. Overlapping but with mixed zygosity"""
  c1 = [v10, v11, v12, v13] = [nv(1, 4, 'CAA', 'C', HOMOZYGOUS),
                               nv(13, 16, 'CAA', 'C', HET_10),
                               nv(20, 23, 'CAA', 'C', HET_10),
                               nv(26, 29, 'CAA', 'C', HOMOZYGOUS)]
  c2 = [v20, v21, v22, v23] = [nv(5, 8, 'CTT', 'C', HOMOZYGOUS),  # Will collide out
                               nv(13, 16, 'CTT', 'C', HET_01),  # Will pass
                               nv(20, 23, 'CTT', 'C', HOMOZYGOUS),  # Will collide out
                               nv(26, 29, 'CTT', 'C', HET_01)]  # Will collide out
  c3 = vr.merge_variants(c1, c2)
  # assert_sequence_equal(c2, [new_variation(1, 4, 'CAA', 'C', HOMOZYGOUS),
  #                   new_variation(13, 16, 'CTT', 'C', HET_10),
  #                   new_variation(13, 16, 'CAA', 'C', HET_01),
  #                   new_variation(20, 23, 'CAA', 'C', HET_10),
  #                   new_variation(26, 29, 'CGG', 'C', HOMOZYGOUS)])
  assert_sequence_equal([v10, v11, v21, v12, v13], c3)


def merge_test7():
  """Skip bug, discovered during simulations"""
  cv1 = nv(3, 7, 'ACTG', 'A', HOMOZYGOUS)
  cv2 = nv(10, 11, 'C', 'A', HOMOZYGOUS)
  c1 = [cv1, cv2]

  dv1 = nv(1, 5, 'ACTG', 'A', HOMOZYGOUS)
  dv2 = nv(7, 8, 'C', 'T', HOMOZYGOUS)  # Collides with previous, should be discarded
  dv3 = nv(15, 16, 'A', 'T', HOMOZYGOUS)  # Too close together, second one should collide out
  dv4 = nv(16, 16, 'C', 'T', HOMOZYGOUS)
  dnv = [dv1, dv2, dv3, dv4]

  c2 = vr.merge_variants(c1, dnv)
  assert_sequence_equal(c2, [cv1, cv2, dv3])


def merge_test8():
  """Tricky variant merge case that led to out of order result"""
  cv1 = nv(26, 27, 'A', 'ACA', HET_01)
  c1 = [cv1]

  dv1 = nv(24, 25, 'A', 'C', HET_01)
  dv2 = nv(25, 26, 'C', 'A', HET_10)
  dnv = [dv1, dv2]

  c2 = vr.merge_variants(c1, dnv)
  assert_sequence_equal(c2, [dv2, cv1], c2)


def copy_missing_chromosomes_test():
  """Copy missing chromosomes."""
  c1 = [v10, v11] = [nv(1, 2, 'A', 'T', HOMOZYGOUS),
                     nv(4, 5, 'C', 'A', HOMOZYGOUS)]
  g1 = {1: c1}
  c2 = [v20, c21] = [nv(3, 4, 'A', 'T', HOMOZYGOUS),
                     nv(6, 7, 'C', 'A', HOMOZYGOUS)]
  g2 = {2: c2}
  vr.copy_missing_chromosomes(g1, g2)
  assert_sequence_equal(g1[1], c1)
  assert_sequence_equal(g1[2], c2)


def merge_genomes_test():
  """Merge genomes."""
  c1 = [v10, v11] = [nv(1, 2, 'A', 'T', HOMOZYGOUS),
                     nv(4, 5, 'C', 'A', HOMOZYGOUS)]
  g1 = {1: c1}
  c2 = [v20, c21] = [nv(3, 4, 'A', 'T', HOMOZYGOUS),
                     nv(6, 7, 'C', 'A', HOMOZYGOUS)]
  g2 = {1: c1, 2: c2}
  g3 = vr.merge_genomes(g1, g2)
  assert_sequence_equal(g3[1], c1)
  assert_sequence_equal(g3[2], c2)
