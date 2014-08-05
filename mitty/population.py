"""
This module implements a haploid genome as a diff with respect to a reference (which never needs to be explicitly set).
The contents of the class, therefore, correspond directly to a very strict version of the VCF and each genome can be
written out as a VCF file.

The internal structure of an attached pair of chromosomes is a list of tuples of the form:

(start, stop, REF, ALT, het)

The module implements the following useful population functions and their supporting functions

crossover(g1, rng, params) -> g2  simulating cross over between copies of chromosomes
fertilization(g1, g2, rng)   -> g3  simulating creation of a child from parents
denovo(g1, ref, rng, params)  -> g3  the only function that requires the original seq - generates denovo mutations

The module implements and uses the following utility functions

parse_vcf(vcf_reader)  -> g1 read in a VCF file and convert it to out haploid genome format
vcf2chrom(vcf_reader)  -> c1 read in one chromosome from a VCF file
write_to_vcf(fname, g1) write the data to a compressed, indexed VCF file

Internally, the only operation that violates sorting order is denovo, so this uses a blist.sortedlist. For everything
else we use a plain Python list as we only append items (and that is O(1) for a list)




get_rngs(seed)  ->  rng  return a list of random number generators that are used by the different functions

The internal genome format is a simple list of tuples almost corresponding to the VCF

(start, stop, REF, ALT)

The following worker functions operate on single copies of chromosomes

merge(c1, c2) -> c3  merge the descriptions of the two chromosomes
                     This figures out any het mutations. The output is sorted and is of the VCF format
                    (POS, POS1, REF, ALT, HET) We only need this when saving back to VCF



Though we could have written a class called genome, I prefer to write in as functional a style as possible for better
code quality.
"""
import numpy
from variation import Variation, HOMOZYGOUS, HET1, HET2


def vcf2chrom(vcf_rdr):
  """Given a vcf reader corresponding to one chromosome, read in the variant descriptions into our format. The result is
  sorted if the vcf file is sorted.
  """
  chrom = []

  for variant in vcf_rdr:
    alt = variant.ALT[0].sequence if variant.ALT[0] is not None else ''
    ref = variant.REF or ''
    start = variant.POS  # Note, we are in VCF coordinates!
    stop = variant.POS + len(ref)
    het = HOMOZYGOUS

    try:
      if variant.samples[0].gt_nums[0] == '0':
        het = HET2
      if variant.samples[0].gt_nums[2] == '0':
        if het == HET2:  # 0/0 means this does not exist in this sample
          continue
        else:
          het = HET1
    except IndexError:  # No genotype info, will assume homozygous
        pass

    chrom += [Variation(start, stop, ref, alt, het)]

  return chrom


def parse_vcf(vcf_rdr, chrom_list):
  """Given a vcf reader load in all the chromosomes."""
  g1 = {}
  for chrom in chrom_list:
    try:
      g1[chrom] = vcf2chrom(vcf_rdr.fetch(chrom, start=0))
    except KeyError:
      g1[chrom] = []

  return g1


def chrom_crossover(c1, crossover_idx):
  """cross_over_idx is a list the same size as c1 with a 1 (indicating a crossover should be done) or 0 (no crossover).
  """
  c2 = []
  for c, idx in zip(c1, crossover_idx):
    if idx == 0:
      c2 += [c]
    elif c.het == HOMOZYGOUS:
      c2 += [c]
    elif c.het == HET1:
      c2 += [c._replace(het=HET2)]
    else:
      c2 += [c._replace(het=HET1)]
  return c2


def crossover_event(g1, crossover_idx):
  return {chrom: chrom_crossover(g1[chrom], crossover_idx[chrom]) for chrom in g1}


def pair_one_chrom(c1, c2, which_copy):
  """Has a resemblance to merge sort BUT WITH GENOMIC VARIANTS!
  Parameters
  ----------
  c1, c2     : tuple list
               (POS, POS2, REF, ALT, HET)
  which_copy : tuple
               e.g. (0, 1) telling us which chromosome of each individual to take

  Returns
  -------
  c3         : tuple list
               (POS, POS2, REF, ALT, HET)
  """
  c1_iter, c2_iter = c1.__iter__(), c2.__iter__()
  c3 = []
  # Try the zipper
  l1, l2 = next(c1_iter, None), next(c2_iter, None)
  while l1 is not None and l2 is not None:

    if (l1.het == HET1 and which_copy[0] == 1) or (l1.het == HET2 and which_copy[0] == 0):
      l1 = next(c1_iter, None)
      continue

    if (l2.het == HET1 and which_copy[1] == 1) or (l2.het == HET2 and which_copy[1] == 0):
      l2 = next(c2_iter, None)
      continue

    if l1._replace(het=HOMOZYGOUS) == l2._replace(het=HOMOZYGOUS):  # Homozygous
      c3 += [l1._replace(het=HOMOZYGOUS)]
      l1, l2 = next(c1_iter, None), next(c2_iter, None)
      continue

    if l1.start <= l2.start:
      c3 += [l1._replace(het=HET1)]
      l1 = next(c1_iter, None)
    else:
      c3 += [l2._replace(het=HET2)]
      l2 = next(c2_iter, None)

  # Now pick up any slack
  while l1 is not None:
    if (l1.het == HET1 and which_copy[0] == 1) or (l1.het == HET2 and which_copy[0] == 0):
      pass
    else:
      c3 += [l1._replace(het=HET1)]
    l1 = next(c1_iter, None)

  while l2 is not None:
    if (l2.het == HET1 and which_copy[1] == 1) or (l2.het == HET2 and which_copy[1] == 0):
      pass
    else:
      c3 += [l2._replace(het=HET2)]
    l2 = next(c2_iter, None)

  return c3


def fertilize_one(g1, g2, which_copy):
  #return {chrom: pair_one_chrom(c1, c2, wc) for (chrom, c1), (_, c2), wc in zip(g1.iteritems(), g2.iteritems(), which_copy)}
  return {chrom: pair_one_chrom(g1[chrom], g2[chrom], which_copy[chrom]) for chrom in g1}


def place_crossovers_on_chrom(c1, hot_spots, rng):
  """Stock cross over generator. Place segments based on gaussian distribution around hot_spots
  The model is a mixture-gaussian model: the hotspot has an amplitude and a width (sd of the gaussian). The amplitude of
  the gaussian at a certain locus determines the probability of the variants there crossing over. The probability of
  crossover at any point is the sum of all Gaussians clipped at +1.0

  Parameters
  ----------
  c1               : list of tuples
                     [(start, stop, <any other payload>) ...]. Should be sorted
  hot_spots        : numpy array
                     n x 3
                     [(center, max_p, width) ... ] 0 <= max_p <= 1
  rng              : random number generator

  Returns
  -------
  idx              : numpy array
                     As used by chrom_crossover
  """
  if hot_spots.shape[0] == 0:
    return [0] * len(c1)  # No hotspots, no crossovers
  #x = numpy.array(c1, dtype=[('st', float), ('a', 'c'), ('b', 'c'), ('c', 'c'), ('d', 'c')])['st']
  x = [v.start for v in c1]
  X, C = numpy.meshgrid(x, hot_spots[:, 0])
  _, A = numpy.meshgrid(x, hot_spots[:, 1])
  _, W = numpy.meshgrid(x, hot_spots[:, 2])
  p = (A * numpy.exp(-(((X - C) / W) ** 2))).sum(axis=0).clip(0, 1.0)
  #return numpy.nonzero(rng.uniform(size=p.size) < p)[0]
  return rng.uniform(size=p.size) < p


def place_crossovers(g1, hot_spots, rng):
  return {chrom: place_crossovers_on_chrom(g1[chrom], hot_spots[chrom], rng) for chrom in g1}


def which_copies(g1, rng):
  """Simple 50/50 coin toss to decide which chromosome copy gets picked during fertilization."""
  return {chrom: wc for chrom, wc in zip(g1, rng.randint(0, 2, size=(len(g1), 2)))}


def get_rngs(seed):
  """There are two rngs needed, one for the cross over hot spots and the other to decide which copy of a chromosome
   goes to a gamete."""
  return [numpy.random.RandomState(seed=sub_seed) for sub_seed in numpy.random.RandomState(seed=seed).randint(100000000, size=2)]


def spawn(g1, g2, hot_spots={}, rngs=[], num_children=2):
  if g1.keys() != g2.keys():
    raise RuntimeError('Two genomes have unequal chromosomes')
  children = []
  for _ in range(num_children):
    g1_cross = crossover_event(g1, place_crossovers(g1, hot_spots, rngs[0]))
    g2_cross = crossover_event(g2, place_crossovers(g2, hot_spots, rngs[0]))
    children.append(fertilize_one(g1_cross, g2_cross, which_copies(g1, rngs[1])))
  return children


def de_novo_population():
  """This uses routines in mutate to introduce denovo mutations ."""