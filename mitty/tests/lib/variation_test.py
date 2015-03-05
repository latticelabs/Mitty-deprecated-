import mitty.lib.variation as vr
from nose.tools import assert_sequence_equal
nvr = vr.Variant
avm = vr.add_novel_variant_to_master
HOM = vr.HOM
HET_01 = vr.HET_01
HET_10 = vr.HET_10


# Some utility functions to make our life easier
def nsv(pos, stop, ref, alt, gt):
  """New sample variant"""
  return vr.SampleVariant(gt, vr.Variant(pos, stop, ref, alt))


def nsp(svs):
  """Given a list of SampleVariants, turn them into a Sample instance"""
  s = vr.Sample('test')
  for sv in svs:
    s.append(sv)
  return s


def avms(svs, ml):
  for sv in svs:
    avm(sv.data, ml)


def addition_test():
  """Add variant to master list"""
  l = {}
  v1 = nvr(1, 2, 'A', 'C')
  assert vr.add_novel_variant_to_master(v1, l) == v1
  assert l[65536] == v1  # Test if index stuff is correct


def addition_test2():
  """Duplicate variant in master list"""
  l = {}
  v1 = nvr(1, 2, 'A', 'C')
  v2 = nvr(1, 2, 'A', 'C')

  assert vr.add_novel_variant_to_master(v1, l) == v1
  assert vr.add_novel_variant_to_master(v2, l) == v1  # Duplicate should be replaced with existing


def addition_test3():
  """Multiple allele detection on addition to master list"""
  l = {}
  v1 = nvr(1, 2, 'A', 'C')
  v2 = nvr(1, 2, 'A', 'G')
  v3 = nvr(1, 2, 'C', 'G')

  assert vr.add_novel_variant_to_master(v1, l) == v1
  assert vr.add_novel_variant_to_master(v2, l) == v2
  assert vr.add_novel_variant_to_master(v3, l) == v3
  assert l[65536] == v1
  assert l[65537] == v2
  assert l[65538] == v3


def sample_test():
  """Appending to Sample linked list"""
  s = vr.Sample('test')
  assert s.label == 'test'
  assert s.head == s.tail
  assert s.length == 0

  sv1 = nsv(1, 4, 'CAA', 'C', HOM)
  sv2 = nsv(10, 13, 'CAA', 'C', HOM)

  s.append(sv1)
  assert s.length == 1
  assert s.head.next == sv1, (s.head, sv1)
  assert s.tail == sv1, (s.tail, sv1)

  s.append(sv2)
  assert s.length == 2
  assert s.head.next == sv1, (s.head, sv1)
  assert s.tail == sv2, (s.tail, sv2)


def sample_test1():
  """Advancing through Sample"""
  s = vr.Sample('test')
  sv1 = nsv(1, 4, 'CAA', 'C', HOM)
  sv2 = nsv(10, 13, 'CAA', 'C', HOM)
  s.append(sv1)
  s.append(sv2)

  assert s.advance() == sv1
  assert s.advance() == sv2
  assert s.advance() is None
  assert s.advance() is None

  s.rewind_cursor()
  assert s.advance() == sv1


def sample_test2():
  """Inserting into Sample linked list"""
  s = vr.Sample('test')
  sv1 = nsv(1, 4, 'CAA', 'C', HOM)
  sv2 = nsv(10, 13, 'CAA', 'C', HOM)
  sv3 = nsv(2, 5, 'CAA', 'C', HOM)
  sv4 = nsv(6, 9, 'CAA', 'C', HOM)

  s.append(sv1)
  s.append(sv2)

  s.advance()
  s.insert(sv3)
  assert_sequence_equal([sv for sv in s], [sv3, sv1, sv2])
  s.insert(sv4)
  assert_sequence_equal([sv for sv in s], [sv3, sv4, sv1, sv2])

  s = vr.Sample('test')
  s.append(sv1)
  s.append(sv2)
  s.advance()
  s.advance()
  s.insert(sv3)
  assert_sequence_equal([sv for sv in s], [sv1, sv3, sv2])

  s = vr.Sample('test')
  s.append(sv1)
  s.append(sv2)
  s.advance()
  s.advance()
  s.advance()
  s.insert(sv3)
  assert_sequence_equal([sv for sv in s], [sv1, sv2, sv3])


def merge_test1():
  """Merge variants, non overlapping, existing first (ED)."""
  l = {}
  v10 = nsv(1, 4, 'CAA', 'C', HOM)
  v20 = nsv(10, 13, 'CAA', 'C', HOM)
  avm(v10.data, l)

  c1 = nsp([v10])
  dnv = iter([v20])
  vr.add_denovo_variants_to_sample(c1, dnv, l)
  assert_sequence_equal(c1.to_list(), [v10, v20])


def merge_test2():
  """Merge variants, non overlapping, denovo first (DE)."""
  l = {}
  v10 = nsv(10, 13, 'CAA', 'C', HOM)
  v20 = nsv(1, 4, 'CAA', 'C', HOM)
  avm(v10.data, l)

  c1 = nsp([v10])
  dnv = iter([v20])
  vr.add_denovo_variants_to_sample(c1, dnv, l)
  assert_sequence_equal(c1.to_list(), [v20, v10])


def merge_test3():
  """Merge variants, non overlapping (EDED)"""
  l = {}
  v10 = nsv(1, 4, 'CAA', 'C', HOM)
  v11 = nsv(13, 16, 'CTT', 'C', HOM)
  v20 = nsv(8, 11, 'CCC', 'C', HOM)
  v21 = nsv(20, 23, 'CAA', 'C', HOM)
  avms([v10, v11], l)

  c1 = nsp([v10, v11])
  dnv = iter([v20, v21])
  vr.add_denovo_variants_to_sample(c1, dnv, l)
  assert_sequence_equal(c1.to_list(), [v10, v20, v11, v21])


def merge_test4a():
  """Merge variants, full overlapping (E-D-)"""
  l = {}
  v10 = nsv(2, 5, 'CAA', 'C', HOM)
  v20 = nsv(2, 5, 'CAA', 'T', HOM)
  avms([v10], l)

  c1 = nsp([v10])
  dnv = iter([v20])
  vr.add_denovo_variants_to_sample(c1, dnv, l)
  assert_sequence_equal(c1.to_list(), [v10])


def merge_test4b():
  """Merge variants, full overlapping, SNP (D-E-)."""
  l = {}
  v10 = nsv(2, 3, 'C', 'G', HET_10)
  v20 = nsv(2, 3, 'C', 'T', HET_10)
  avms([v10], l)

  c1 = nsp([v10])
  dnv = iter([v20])
  vr.add_denovo_variants_to_sample(c1, dnv, l)
  assert_sequence_equal(c1.to_list(), [v10])


# Important test - killed a nasty logic bug
def merge_test4c():
  """Merge variants, full overlapping, with non-colliding preceder"""
  l = {}
  v10 = nsv(1, 2, 'G', 'C', HET_01)
  v11 = nsv(2, 3, 'A', 'C', HET_10)
  v20 = nsv(2, 3, 'A', 'T', HET_10)
  avms([v10, v11], l)

  c1 = nsp([v10, v11])
  dnv = iter([v20])
  vr.add_denovo_variants_to_sample(c1, dnv, l)
  assert_sequence_equal(c1.to_list(), [v10, v11])


def merge_test4():
  """Merge variants, overlapping (E-D). D will collide and will be rejected"""
  l = {}
  v10 = nsv(1, 4, 'CAA', 'C', HOM)
  v20 = nsv(2, 5, 'CCC', 'C', HOM)
  avms([v10], l)

  c1 = nsp([v10])
  dnv = iter([v20])
  vr.add_denovo_variants_to_sample(c1, dnv, l)
  assert_sequence_equal(c1.to_list(), [v10])


def merge_test5():
  """Merge variants, overlapping (D-E). D will collide and will be rejected"""
  l = {}
  v10 = nsv(2, 5, 'CCC', 'C', HOM)
  v20 = nsv(1, 4, 'CAA', 'C', HOM)
  avms([v10], l)

  c1 = nsp([v10])
  dnv = iter([v20])
  vr.add_denovo_variants_to_sample(c1, dnv, l)
  assert_sequence_equal(c1.to_list(), [v10])


def merge_test6():
  """Merge variants, overlapping heterozygous. Overlapping but with mixed zygosity"""
  l = {}
  v10 = nsv(1, 4, 'CAA', 'C', HOM)
  v11 = nsv(13, 16, 'CAA', 'C', HET_10)
  v12 = nsv(20, 23, 'CAA', 'C', HET_10)
  v13 = nsv(26, 29, 'CAA', 'C', HOM)
  avms([v10, v11, v12, v13], l)

  v20 = nsv(5, 8, 'CTT', 'C', HOM)  # Will collide out
  v21 = nsv(13, 16, 'CTT', 'C', HET_01)
  v22 = nsv(20, 23, 'CTT', 'C', HOM)
  v23 = nsv(26, 29, 'CTT', 'C', HET_01)

  c1 = nsp([v10, v11, v12, v13])
  dnv = iter([v20, v21, v22, v23])
  vr.add_denovo_variants_to_sample(c1, dnv, l)
  assert_sequence_equal(c1.to_list(), [v10, v11, v21, v12, v13])


def merge_test7():
  """Skip bug, discovered during simulations"""
  l = {}
  v10 = nsv(3, 7, 'ACTG', 'A', HOM)
  v11 = nsv(10, 11, 'C', 'A', HOM)
  avms([v10, v11], l)

  v20 = nsv(1, 5, 'ACTG', 'A', HOM)
  v21 = nsv(7, 8, 'C', 'T', HOM)    # Collides with previous, should be discarded
  v22 = nsv(15, 16, 'A', 'T', HOM)  # Too close together, second one should collide out
  v23 = nsv(16, 16, 'C', 'T', HOM)

  c1 = nsp([v10, v11])
  dnv = iter([v20, v21, v22, v23])
  vr.add_denovo_variants_to_sample(c1, dnv, l)
  assert_sequence_equal(c1.to_list(), [v10, v11, v22])


def merge_test8():
  """Tricky data merge case that led to out of order result"""
  l = {}
  v10 = nsv(26, 27, 'A', 'ACA', HET_01)
  avms([v10], l)

  v20 = nsv(24, 25, 'A', 'C', HET_01)
  v21 = nsv(25, 26, 'C', 'A', HET_10)

  c1 = nsp([v10])
  dnv = iter([v20, v21])
  vr.add_denovo_variants_to_sample(c1, dnv, l)
  assert_sequence_equal(c1.to_list(), [v21, v10])


# def copy_missing_chromosomes_test():
#   """Copy missing chromosomes."""
#   l = {}
#   c1 = [v10, v11] = [cgt(l, 1, 2, 'A', 'T', HOMOZYGOUS),
#                      cgt(l, 4, 5, 'C', 'A', HOMOZYGOUS)]
#   g1 = {1: c1}
#   c2 = [v20, c21] = [cgt(l, 3, 4, 'A', 'T', HOMOZYGOUS),
#                      cgt(l, 6, 7, 'C', 'A', HOMOZYGOUS)]
#   g2 = {2: c2}
#   vr.copy_missing_chromosomes(g1, g2)
#   assert_sequence_equal(g1[1], c1)
#   assert_sequence_equal(g1[2], c2)
#
#
# def merge_genomes_test():
#   """Merge genomes."""
#   l = {}
#   c1 = [v10, v11] = [cgt(l, 1, 2, 'A', 'T', HOMOZYGOUS),
#                      cgt(l, 4, 5, 'C', 'A', HOMOZYGOUS)]
#   g1 = {1: c1}
#   c2 = [v20, c21] = [cgt(l, 3, 4, 'A', 'T', HOMOZYGOUS),
#                      cgt(l, 6, 7, 'C', 'A', HOMOZYGOUS)]
#   g2 = {1: c1, 2: c2}
#   g3 = vr.merge_genomes(g1, g2, l)
#   assert_sequence_equal(g3[1], c1)
#   assert_sequence_equal(g3[2], c2)
