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
from variation import vcopy, HOMOZYGOUS, HET1, HET2
import mitty.denovo


def chrom_crossover(c1, crossover_idx):
  """cross_over_idx is a list the same size as c1 with a 1 (indicating a crossover should be done) or 0 (no crossover).
  """
  c2 = []
  for c, idx in zip(c1, crossover_idx):
    new_c = vcopy(c)  # Valid for no crossover or HOMOZYGOUS
    if idx == 1 and c.het != HOMOZYGOUS:
      if c.het == HET1:
        new_c.het = HET2
      else:
        new_c.het = HET1
    c2 += [new_c]
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

    if vcopy(l1, het=HOMOZYGOUS) == vcopy(l2, het=HOMOZYGOUS):  # Homozygous
      c3 += [vcopy(l1, het=HOMOZYGOUS)]
      l1, l2 = next(c1_iter, None), next(c2_iter, None)
      continue

    if l1.POS <= l2.POS:
      c3 += [vcopy(l1, het=HET1)]
      l1 = next(c1_iter, None)
    else:
      c3 += [vcopy(l2, het=HET2)]
      l2 = next(c2_iter, None)

  # Now pick up any slack
  while l1 is not None:
    if (l1.het == HOMOZYGOUS) or (l1.het == HET1 and which_copy[0] == 0) or (l1.het == HET2 and which_copy[0] == 1):
      c3 += [vcopy(l1, het=HET1)]
    l1 = next(c1_iter, None)

  while l2 is not None:
    if (l2.het == HOMOZYGOUS) or (l2.het == HET1 and which_copy[1] == 0) or (l2.het == HET2 and which_copy[1] == 1):
      c3 += [vcopy(l2, het=HET2)]
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
  x = [v.POS for v in c1]
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
  """Creates independent RNGs for our random processes from the master seed given."""
  rng_names = ['cross_over', 'chrom_copy', 'couple_chose']
  return {k: numpy.random.RandomState(seed=sub_seed)
          for k,sub_seed in zip(rng_names, numpy.random.RandomState(seed=seed).randint(100000000, size=len(rng_names)))}


def spawn(g1, g2, hot_spots={}, rngs={}, num_children=2, ref_fp=None, models=[]):
  """If you don't want denovo mutations don't send in ref_fp and models"""
  if g1.keys() != g2.keys():
    raise RuntimeError('Two genomes have unequal chromosomes')
  children = []
  for _ in range(num_children):
    g1_cross = crossover_event(g1, place_crossovers(g1, hot_spots, rngs['cross_over']))
    g2_cross = crossover_event(g2, place_crossovers(g2, hot_spots, rngs['cross_over']))
    if ref_fp is not None and len(models) > 0:  # We want denovo mutations
      mitty.denovo.add_multiple_variant_models_to_genome(ref_fp, models, g1_cross)
      mitty.denovo.add_multiple_variant_models_to_genome(ref_fp, models, g2_cross)
    children.append(fertilize_one(g1_cross, g2_cross, which_copies(g1, rngs['chrom_copy'])))
  return children


def de_novo_population(ref_fp, models=[], size=10):
  """This uses variant plugins and denovo.py functions to generate a population of highly differentiated individuals"""
  return [mitty.denovo.create_denovo_genome(ref_fp, models) for _ in range(size)]


def one_generation(pop, hot_spots={}, rngs={}, num_children_per_couple=2, ref_fp=None, models=[]):
  """We take pairs from the population without replacement, cross over the chromosomes, sprinkle in denovo mutations
  then shuffle the chromosomes during fertilization. If the population size is odd, we discard one parent"""
  assert len(pop) >= 2, 'Need at least two parents'
  mates = rngs['couple_chose'].permutation(len(pop))  # This needs to be an even number
  children = []
  parents = []
  for n in range(len(mates))[::2]:
    children += spawn(pop[mates[n]], pop[mates[n + 1]], hot_spots=hot_spots, rngs=rngs, num_children=num_children_per_couple,
                      ref_fp=ref_fp, models=models)
    parents += [(mates[n], mates[n + 1])]
  return children, parents


def population_simulation(ref_fp, denovo_models=[], initial_size=10,
                          hot_spots={}, rngs={}, num_children_per_couple=2, ss_models=[],
                          num_generations=10,
                          store_all_generations=True):
  # The initial population
  pop = de_novo_population(ref_fp, denovo_models, size=initial_size)

  generation = [pop]
  parent_list = []
  # The generational loop
  this_gen = pop
  for n in range(num_generations):
    children, parents = one_generation(this_gen, hot_spots=hot_spots, rngs=rngs,
                                       num_children_per_couple=num_children_per_couple,
                                       ref_fp=ref_fp, models=ss_models)
    parent_list.append(parents)
    if store_all_generations:
      generation.append(children)
    this_gen = children
  return generation, parent_list

#def save_population(pop, gen):

