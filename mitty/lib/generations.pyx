"""This is to simulate sexual reproduction in a population of individuals. This module takes care of the higher level
methods and is designed to be independent on Variation. The idea is that plugins to determine fitness, place crossovers
etc, do the actual crossing-over, will operate on the genome side (which is carried as a pointer in the Sample structure)
but we don't need any of the details to do the overall population simulation.

1. A founder population of completely unrelated samples is created
2. A mate picker algorithm picks pairs of parents from the population for mating
3. A breeder algorithm creates children for each parent
4. A culler algorithm removes unfit children and keeps the population within limits
5. A fitness algorithm is used by several of the other algorithms to give a fitness value to a genome
6. A fracture algorithm can be applied to split the population


The stock plugins are described later below

1. mating plugin: randomly picks pairs of mates from the pool. A picked pair is rejected if they have the same parent
N generations back. The number of children is determined by the combined fitness values of the parents

2. cross-over: has a flat probability of cross over points. If you find hotspots, your algorithm is flawed

3. fitness: has an (arbitrary) map of regions where mutations are deleterious and where they are beneficial and where
mutations are recessive and dominant

4. fracture: break the given population into two separate parts
"""
import random

import mitty.lib.variation as vr  # Only needed for the genomify function. Everything else is independent
import mitty.lib.util as mutil


cdef class Sample:
  """This represents one genome sample.

  Attributes:
    generation - which generation is this sample
    serial     - what number in this generation
    name       - string for the form 'generation:serial'
    parents    - list of references to parent Samples
    genome     - the genome dictionary of this sample
    fitness    - computed fitness value of the genome. [-1.0,1.0]
  """
  cdef public:
    str name
    int generation, serial, hash
    float fitness
    Sample p1, p2
    dict genome


  def __cinit__(self, generation, serial):
    self.name = '{:d}:{:d}'.format(generation, serial)  # Our unique sample name
    self.generation, self.serial, self.p1, self.p2 = generation, serial, None, None
    self.genome, self.fitness = {}, 0.0  # Fitness is filled by an external function.
    self.hash = self.generation ^ self.serial

  def __hash__(self):
    return self.hash

  def __richcmp__(self, Sample other, int op):
    if op == 2:
      return (self.serial == other.serial) and (self.generation == other.generation)
    elif op == 3:
      return not ((self.serial == other.serial) and (self.generation == other.generation))

  def __repr__(self):
    return '{:s} ({:1.3f}) [{:s}, {:s}]'.format(self.name, self.fitness, self.p1.name, self.p2.name)


def create_initial_generation(generation_size):
  """:returns: List of Sample with generation set to 0"""
  return [Sample(0, n) for n in range(generation_size)]


def create_next_generation(parent_generation, fitness, mater, breeder, culler):
  """Wrapper round plugins to select mates from the population and create children.

  :param parent_generation: list of Sample representing parent population
  :param mater: mate picker function
  :param breeder: child generator function
  :fitness: fitness class
  :returns: list of lists of children. Children are grouped by parent pair, though each child has their parent filled out
  """
  gen = parent_generation[0].generation if len(parent_generation) else 0
  mate_list = mater.pick_mates(parent_generation)
  children = breeder.breed(mate_list, gen + 1)
  fitness.fitness(children)
  culler.cull(children)
  return children


cpdef genomify(list children, cross_model):
  """For each child, find the parents and pair their genomes"""
  cdef Sample ch
  for ch in children:
    ch.genome = vr.pair_samples(ch.p1.genome, cross_model.cross_overs(ch.p1.genome), cross_model.copies(ch.p1.genome),
                                ch.p2.genome, cross_model.cross_overs(ch.p2.genome), cross_model.copies(ch.p2.genome))



# ------------------ STOCK POPULATION SIM PLUGINS ----------------------------- #


class StockFitness:
  """Simple fitness computation. Rewards variations in second quarter of genome and penalized variations in third q"""
  def __init__(self, ref):
    """:param ref: list of tuples as returned by io.load_multi_fasta_gz"""
    l = self.chromosome_lengths = [seq[1] for seq in ref]
    self.l0 = [0.25 * ll for ll in l]
    self.l1 = [0.5 * ll for ll in l]
    self.l2 = [0.75 * ll for ll in l]

  def sample_fitness(self, sample):
    """Modifies fitness in place"""
    g = sample.genome
    fit = 0
    v_count = 0
    for ch, chrom in g:
      for gtv in chrom:
        if self.l0 < gtv.pos < self.l1:
          fit += 1.0
        elif self.l1 < gtv.pos < self.l2:
          fit -= 1.0
        else:
          continue
        v_count += 1
    sample.fitness = fit / v_count if v_count else 0

  def fitness(self, generation):
    for g in generation:
      self.sample_fitness(g)
    #map(self.sample_fitness, generation)


class StockMater:
  """Simple stock mating algorithm. Shuffles parent list and returns mate pairs that are not related to N generations."""
  def __init__(self, incest_generations=2, master_seed=2):
    """
    :param incest_generations: number of generations to go back for relatedness test
    :param master_seed: Integer seed for RNG"""
    self.incest_generations, self.rnd = incest_generations, random.Random(master_seed)

  def incestuous_mating(self, s1, s2):
    """Return False if the samples are not related up to the generations indicated. (Find if there is an LCA within x levels)
    :param s1, s2: samples to test
    :returns boolean
    """
    parent_list = [p for s in [s1, s2] for p in [s.p1, s.p2] if p is not None]
    for n in range(self.incest_generations):
      if len(set(parent_list)) < len(parent_list):
        return True
      parent_list = [p for s in parent_list for p in [s.p1, s.p2]]
    return False

  def mate_in_sequence(self, parent_list):
    """Simple mating algorithm that pops pairs of parents from the list, tests for incest, and then adds legal pairs to a
    list and illegal pairs to another list
    :param parent_list: list of Samples that will be parents. Items will be popped from this list, just so you know
    :returns mating_list - list of pairs of legal mates
             incest_list - two lists of samples. Corresponding samples in the lists are related
    """
    mating_list = []
    incest_list = [[], []]  # The folks we couldn't mate, kept in separate lists to increase their chances of mating
    while len(parent_list) > 1:
      s1 = parent_list.pop()
      s2 = parent_list.pop()
      if self.incestuous_mating(s1, s2):
        incest_list[0] += [s1]
        incest_list[1] += [s2]
        continue
      mating_list += [[s1, s2]]
    return mating_list, incest_list

  def pick_mates(self, parent_list):
    """Simple mate picker that shuffles the parent_list and then pops pairs of parents off the shuffled list.
    :param parent_list: list of Samples that will be parents.
    """
    _parent_list = [p for p in parent_list]  # Better not mess with the original
    self.rnd.shuffle(_parent_list)
    mating_list, incest_list = self.mate_in_sequence(_parent_list)

    # Make an attempt on the rejects
    mating_list2, _ = self.mate_in_sequence(incest_list[0])
    mating_list3, _ = self.mate_in_sequence(incest_list[1])

    return mating_list + mating_list2 + mating_list3


class StockBreeder:
  """Stock breeding algorithm, creates children based on average fitness of parents."""
  def __init__(self, child_factor=2.0, master_seed=3):
    """:param child_factor: Average number of children for average fitness parents"""
    self.child_factor, self.rng = child_factor, random.Random(master_seed)

  def breed(self, mate_list, generation):
    """
    :param mate_list: list of pairs of parents from which we generate children
    :param generation: int indicating which generation the children will be
    :returns: list of Samples"""
    children = []
    cnt = 0
    for p1, p2 in mate_list:
      avg_fit = (p1.fitness + p2.fitness) / 2.0
      n_children = int(max(0, self.rng.gauss(self.child_factor * (avg_fit + 1.0), 1)) + .5)
      these_children = []
      for n in range(n_children):
        this_child = Sample(generation, cnt)
        this_child.p1, this_child.p2 = p1, p2
        these_children += [this_child]
        cnt += 1
      children += these_children
    return children


class StockCuller:
  """Given a maximum population size, cut down the least fit children"""
  def __init__(self, max_pop_size):
    self.max_pop_size = max_pop_size

  def cull(self, population):
    while len(population) > self.max_pop_size:
      population.remove(min(population, key=lambda x: x.fitness))


# ------------------ STOCK POPULATION SIM PLUGIN ----------------------------- #


class StockCrossOver:  # This is more closely affiliated with variations than generations
  """Places crossover points at any point with equal probability."""
  def __init__(self, ref, p=1e-8, master_seed=3):
    """
    :param ref: list of tuples as returned by fasta loader
    :param p: probability of crossover occuring at any base
    :param master_seed: seed for rngs
    """
    self.l = self.chromosome_lengths = [seq[1] for seq in ref]
    self.pos_rng, self.copy_rng = mutil.initialize_rngs(master_seed, 2)
    self.p = p

  def cross_overs(self, g):
    """
    :param g: dict of Chromosome objects in case we want to use the actual sequence data to place cross overs
    :returns: dict of crossover positions as needed by pair_chromosomes
    """
    return {n+1: mutil.place_poisson(self.pos_rng, self.p, self.l[n]).tolist() for n in range(len(self.l))}

  def copies(self, g):
    """
    :param g: dict of Chromosome objects in case we want to use the actual sequence data to place cross overs
    :returns: dict of crossover positions as needed by pair_chromosomes
    """
    cpy = self.copy_rng.randint(2, size=len(self.l))
    return {n+1: cpy[n] for n in range(len(self.l))}