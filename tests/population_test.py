import numpy.testing
from mitty.variation import Variation
from mitty.population import *


def chrom_crossover_test():
  c1 = [
      Variation(1, 2, 'C', 'CAA', HET2),
      Variation(3, 6, 'CAG', 'C', HET1),
      Variation(7, 8, 'G', 'T', HET2),
      Variation(9, 12, 'GTT', '', HET1),
      Variation(13, 16, 'GTT', 'TTG', HOMOZYGOUS),
      Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS)
  ]
  crossover_idx = [1, 1, 1, 0, 1, 0]
  correct_chrom = [
      Variation(1, 2, 'C', 'CAA', HET1),
      Variation(3, 6, 'CAG', 'C', HET2),
      Variation(7, 8, 'G', 'T', HET1),
      Variation(9, 12, 'GTT', '', HET1),
      Variation(13, 16, 'GTT', 'TTG', HOMOZYGOUS),
      Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS)
  ]
  c1_c = chrom_crossover(c1, crossover_idx)
  assert correct_chrom == c1_c, c1_c
  assert c1[0] == (1, 2, 'C', 'CAA', HET2)  # We shouldn't be changing our original


def crossover_event_test():
  c1 = [
    Variation(1, 2, 'C', 'CAA', HET2),
    Variation(3, 6, 'CAG', 'C', HET1),
    Variation(7, 8, 'G', 'T', HET2)
  ]
  c2 = [
    Variation(2, 5, 'CAG', 'C', HET1),
    Variation(7, 8, 'G', 'T', HET2)
  ]
  g1 = {1: c1, 2: c2}
  crossover_idx = {1: [1, 0, 1], 2: [0, 1]}

  correct_g = {
    1: [
      Variation(1, 2, 'C', 'CAA', HET1),
      Variation(3, 6, 'CAG', 'C', HET1),
      Variation(7, 8, 'G', 'T', HET1)
    ],
    2: [
      Variation(2, 5, 'CAG', 'C', HET1),
      Variation(7, 8, 'G', 'T', HET1)
    ]
  }
  assert correct_g == crossover_event(g1, crossover_idx)


def pair_one_chrom_test():
  """Merging: test c1 longer than c2"""
  c1 = [
    Variation(1, 2, 'C', 'CAA', HET2),
    Variation(3, 6, 'CAG', 'C', HET1),
    Variation(17, 18, 'G', 'T', HET2),
    Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS),
  ]
  c2 = [
    Variation(7, 8, 'G', 'T', HET2),
  ]
  which_copy = (1, 0)
  correct_pairing = [
    Variation(1, 2, 'C', 'CAA', HET1),
    Variation(17, 18, 'G', 'T', HET1),
    Variation(23, 26, 'GTT', 'TTG', HET1),
  ]
  assert correct_pairing == pair_one_chrom(c1, c2, which_copy), pair_one_chrom(c1, c2, which_copy)


def pair_one_chrom_test2():
  """Merging: test c2 longer than c1"""
  c1 = [
    Variation(1, 2, 'C', 'CAA', HET2)
  ]
  c2 = [
    Variation(7, 8, 'G', 'T', HET1),
    Variation(17, 18, 'G', 'T', HET1),
  ]
  which_copy = (1, 0)
  correct_pairing = [
    Variation(1, 2, 'C', 'CAA', HET1),
    Variation(7, 8, 'G', 'T', HET2),
    Variation(17, 18, 'G', 'T', HET2),
  ]
  assert correct_pairing == pair_one_chrom(c1, c2, which_copy), pair_one_chrom(c1, c2, which_copy)


def pair_one_chrom_test3():
  """Merging: Comprehensive test"""
  c1 = [
    Variation(1, 2, 'C', 'CAA', HET2),  # Tests homozygosity
    Variation(3, 6, 'CAG', 'C', HET1),  # Tests both variants are not on the copies chosen
    Variation(17, 18, 'G', 'T', HET2),  # Test zipper (several variants should come from c2 before we get to this)
    Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS),  # Tests handling of homozygous variants
    Variation(29, 30, 'T', 'G', HOMOZYGOUS)  # Tests unequal var list lengths
  ]
  c2 = [
    Variation(1, 2, 'C', 'CAA', HET1),
    Variation(3, 6, 'CAG', 'C', HET2),
    Variation(7, 8, 'G', 'T', HET2),
    Variation(9, 12, 'GTT', '', HET1),
    Variation(15, 18, 'GTT', 'G', HET1),  # Test partial overlap on different copies (should not affect each other)
    Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS)
  ]
  which_copy = (1, 0)
  correct_pairing = [
    Variation(1, 2, 'C', 'CAA', HOMOZYGOUS),  # Tests homozygosity
    Variation(9, 12, 'GTT', '', HET2),
    Variation(15, 18, 'GTT', 'G', HET2),  # Test partial overlap on different copies (should not affect each other)
    Variation(17, 18, 'G', 'T', HET1),  # Test zipper (several variants should come from c2 before we get to this)
    Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS),  # Tests handling of homozygous variants
    Variation(29, 30, 'T', 'G', HET1)  # Tests unequal var list lengths
  ]
  assert correct_pairing == pair_one_chrom(c1, c2, which_copy)


def fertilize_one_test():
  c11 = [
    Variation(1, 2, 'C', 'CAA', HET1)
  ]
  c12 = [
    Variation(7, 8, 'G', 'T', HET2),
    Variation(17, 18, 'G', 'T', HET2),
  ]
  g1 = {1: c11, 2: c12}

  c21 = [
    Variation(1, 2, 'C', 'CAA', HET2)
  ]
  c22 = [
    Variation(7, 8, 'G', 'T', HET1),
    Variation(17, 18, 'G', 'T', HET1),
  ]
  g2 = {1: c21, 2: c22}
  which_copy = {1: (0, 1), 2: (1, 0)}

  correct_g3 = {
    1: [
      Variation(1, 2, 'C', 'CAA', HOMOZYGOUS)
    ],
    2: [
      Variation(7, 8, 'G', 'T', HOMOZYGOUS),
      Variation(17, 18, 'G', 'T', HOMOZYGOUS),
    ]
  }
  assert correct_g3 == fertilize_one(g1, g2, which_copy)


def place_crossovers_on_chrom_test():
  """Cross over location generator"""
  c1 = [
    Variation(1, 2, 'C', 'CAA', HET2),
    Variation(3, 6, 'CAG', 'C', HET1),
    Variation(17, 18, 'G', 'T', HET2),
    Variation(23, 26, 'GTT', 'TTG', HOMOZYGOUS),
    Variation(29, 30, 'T', 'G', HOMOZYGOUS)
  ]

  hot_spots = numpy.array([])  # No hotspots, no crossover
  numpy.testing.assert_array_equal(numpy.array([0, 0, 0, 0, 0]),  # Hot spot is narrow and over first variant
                                   place_crossovers_on_chrom(c1, hot_spots, numpy.random.RandomState(seed=1)))

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

