import numpy.testing
from mitty.population import *


def vcf2chrom_test():
  import vcf, io

  #
  #  INS     DEL--     SNP      DEL---     INV-----
  #  0    1  2 3 4  5  6    7   8 9 10  11 12 13 14
  #  ------  -----     ------   ------     --------

  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t1\t.\tC\tCAA\t100\tPASS\t.\tGT\t0/1\n"
    "1\t3\t.\tCAG\tC\t100\tPASS\t.\tGT\t1/0\n"
    "1\t7\t.\tG\tT\t100\tPASS\t.\tGT\t0/1\n"
    "1\t9\t.\tGTT\t.\t100\tPASS\t.\tGT\t1/0\n"
    "1\t13\t.\tGTT\tTTG\t100\tPASS\t.\tGT\t1/1\n"
  )

  correct_chrom = [
      (1, 2, 'C', 'CAA', HET2),
      (3, 6, 'CAG', 'C', HET1),
      (7, 8, 'G', 'T', HET2),
      (9, 12, 'GTT', '', HET1),
      (13, 16, 'GTT', 'TTG', HOMOZYGOUS)
    ]

  chrom = vcf2chrom(vcf.Reader(fsock=io.BytesIO(vcf_str)))
  assert chrom == correct_chrom


def chrom_crossover_test():
  c1 = [
      (1, 2, 'C', 'CAA', HET2),
      (3, 6, 'CAG', 'C', HET1),
      (7, 8, 'G', 'T', HET2),
      (9, 12, 'GTT', '', HET1),
      (13, 16, 'GTT', 'TTG', HOMOZYGOUS),
      (23, 26, 'GTT', 'TTG', HOMOZYGOUS)
    ]
  crossover_idx = [1, 1, 1, 0, 1, 0]
  correct_chrom = [
      (1, 2, 'C', 'CAA', HET1),
      (3, 6, 'CAG', 'C', HET2),
      (7, 8, 'G', 'T', HET1),
      (9, 12, 'GTT', '', HET1),
      (13, 16, 'GTT', 'TTG', HOMOZYGOUS),
      (23, 26, 'GTT', 'TTG', HOMOZYGOUS)
    ]
  c1_c = chrom_crossover(c1, crossover_idx)
  assert correct_chrom == c1_c, c1_c
  assert c1[0] == (1, 2, 'C', 'CAA', HET2)  # We shouldn't be changing our original


def crossover_event_test():
  c1 = [
    (1, 2, 'C', 'CAA', HET2),
    (3, 6, 'CAG', 'C', HET1),
    (7, 8, 'G', 'T', HET2)
  ]
  c2 = [
    (2, 5, 'CAG', 'C', HET1),
    (7, 8, 'G', 'T', HET2)
  ]
  g1 = {1: c1, 2: c2}
  crossover_idx = {1: [1, 0, 1], 2: [0, 1]}

  correct_g = {
    1: [
      (1, 2, 'C', 'CAA', HET1),
      (3, 6, 'CAG', 'C', HET1),
      (7, 8, 'G', 'T', HET1)
    ],
    2: [
      (2, 5, 'CAG', 'C', HET1),
      (7, 8, 'G', 'T', HET1)
    ]
  }
  assert correct_g == crossover_event(g1, crossover_idx)


def pair_one_chrom_test():
  """Merging: test c1 longer than c2"""
  c1 = [
    (1, 2, 'C', 'CAA', HET2),
    (3, 6, 'CAG', 'C', HET1),
    (17, 18, 'G', 'T', HET2),
    (23, 26, 'GTT', 'TTG', HOMOZYGOUS),
  ]
  c2 = [
    (7, 8, 'G', 'T', HET2),
  ]
  which_copy = (1, 0)
  correct_pairing = [
    (1, 2, 'C', 'CAA', HET1),
    (17, 18, 'G', 'T', HET1),
    (23, 26, 'GTT', 'TTG', HET1),
  ]
  assert correct_pairing == pair_one_chrom(c1, c2, which_copy), pair_one_chrom(c1, c2, which_copy)


def pair_one_chrom_test2():
  """Merging: test c2 longer than c1"""
  c1 = [
    (1, 2, 'C', 'CAA', HET2)
  ]
  c2 = [
    (7, 8, 'G', 'T', HET1),
    (17, 18, 'G', 'T', HET1),
  ]
  which_copy = (1, 0)
  correct_pairing = [
    (1, 2, 'C', 'CAA', HET1),
    (7, 8, 'G', 'T', HET2),
    (17, 18, 'G', 'T', HET2),
  ]
  assert correct_pairing == pair_one_chrom(c1, c2, which_copy), pair_one_chrom(c1, c2, which_copy)


def pair_one_chrom_test3():
  """Merging: Comprehensive test"""
  c1 = [
    (1, 2, 'C', 'CAA', HET2),  # Tests homozygosity
    (3, 6, 'CAG', 'C', HET1),  # Tests both variants are not on the copies chosen
    (17, 18, 'G', 'T', HET2),  # Test zipper (several variants should come from c2 before we get to this)
    (23, 26, 'GTT', 'TTG', HOMOZYGOUS),  # Tests handling of homozygous variants
    (29, 30, 'T', 'G', HOMOZYGOUS)  # Tests unequal var list lengths
  ]
  c2 = [
    (1, 2, 'C', 'CAA', HET1),
    (3, 6, 'CAG', 'C', HET2),
    (7, 8, 'G', 'T', HET2),
    (9, 12, 'GTT', '', HET1),
    (15, 18, 'GTT', 'G', HET1),  # Test partial overlap on different copies (should not affect each other)
    (23, 26, 'GTT', 'TTG', HOMOZYGOUS)
  ]
  which_copy = (1, 0)
  correct_pairing = [
    (1, 2, 'C', 'CAA', HOMOZYGOUS),  # Tests homozygosity
    (9, 12, 'GTT', '', HET2),
    (15, 18, 'GTT', 'G', HET2),  # Test partial overlap on different copies (should not affect each other)
    (17, 18, 'G', 'T', HET1),  # Test zipper (several variants should come from c2 before we get to this)
    (23, 26, 'GTT', 'TTG', HOMOZYGOUS),  # Tests handling of homozygous variants
    (29, 30, 'T', 'G', HET1)  # Tests unequal var list lengths
  ]
  assert correct_pairing == pair_one_chrom(c1, c2, which_copy)


def fertilize_one_test():
  c11 = [
    (1, 2, 'C', 'CAA', HET1)
  ]
  c12 = [
    (7, 8, 'G', 'T', HET2),
    (17, 18, 'G', 'T', HET2),
  ]
  g1 = {1: c11, 2: c12}

  c21 = [
    (1, 2, 'C', 'CAA', HET2)
  ]
  c22 = [
    (7, 8, 'G', 'T', HET1),
    (17, 18, 'G', 'T', HET1),
  ]
  g2 = {1: c21, 2: c22}
  which_copy = {1: (0, 1), 2: (1, 0)}

  correct_g3 = {
    1: [
      (1, 2, 'C', 'CAA', HOMOZYGOUS)
    ],
    2: [
      (7, 8, 'G', 'T', HOMOZYGOUS),
      (17, 18, 'G', 'T', HOMOZYGOUS),
    ]
  }
  assert correct_g3 == fertilize_one(g1, g2, which_copy)


def place_crossovers_on_chrom_test():
  """Cross over location generator"""
  c1 = [
    (1, 2, 'C', 'CAA', HET2),  # Tests homozygosity
    (3, 6, 'CAG', 'C', HET1),  # Tests both variants are not on the copies chosen
    (17, 18, 'G', 'T', HET2),  # Test zipper (several variants should come from c2 before we get to this)
    (23, 26, 'GTT', 'TTG', HOMOZYGOUS),  # Tests handling of homozygous variants
    (29, 30, 'T', 'G', HOMOZYGOUS)  # Tests unequal var list lengths
  ]

  hot_spots = numpy.array([[1, 1, .5]])
  numpy.testing.assert_array_equal(numpy.array([1, 0, 0, 0, 0]),  # Hot spot is narrow and over first variant
                                   place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))

  hot_spots = numpy.array([[17, 1, .5]])
  numpy.testing.assert_array_equal(numpy.array([0, 0, 1, 0, 0]),  # Hot spot is narrow and over third variant
                                   place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))

  hot_spots = numpy.array([[1, 1, .5], [17, 1, .5]])
  numpy.testing.assert_array_equal(numpy.array([1, 0, 1, 0, 0]),  # Two narrow hot spots over first and third variants
                                   place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))

  hot_spots = numpy.array([[23, 1, 100], [17, 1, 100]])
  numpy.testing.assert_array_equal(numpy.array([1, 1, 1, 1, 1]),  # Super broad hotspots, covers all
                                   place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))

  # A proper test is actually to do something like this (with c1 set as above)
  # import pylab
  # hot_spots = numpy.array([[3, 1, 2], [23, 1, 5]])
  # rng = numpy.random.RandomState(seed=1)
  # data = numpy.concatenate([numpy.nonzero(place_crossovers_on_chrom(c1, hot_spots, rng))[0] for n in range(100)])
  # pylab.hist(data)
  # pylab.plot()
  # The denser your variant structure, the clearer is the sum of gaussian model


def spawn_test():
  pass

