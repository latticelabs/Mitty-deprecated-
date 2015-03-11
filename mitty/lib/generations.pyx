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


class Sample:
  """This represents one genome sample.

  Attributes:
    generation - which generation is this sample
    serial     - what number in this generation
    name       - string for the form 'generation:serial'
    parents    - list of references to parent Samples
    genome     - the genome dictionary of this sample
    fitness    - computed fitness value of the genome. [-1.0,1.0]
  """
  def __init__(self, generation, serial):
    self.name = '{:d}:{:d}'.format(generation, serial)  # Our unique sample name
    self.generation, self.serial, self.parents = generation, serial, []
    self.genome, self.fitness = {}, 0  # Fitness is filled by an external function.

  def __hash__(self):
    return self.name.__hash__()

  def __eq__(self, other):
    return self.__hash__() == other.__hash__()

  def __repr__(self):
    return '{:s} ({:d})'.format(self.name, self.fitness)


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
  children = breeder.breed(mate_list, gen)
  fitness.fitness(children)
  culler.cull(children)
  return children


# ------------------ STOCK POPULATION SIM PLUGINS ----------------------------- #


class StockFitness:
  """Simple fitness computation. Rewards variations in second quarter of genome and penalized variations in third q"""
  def __init__(self, ref):
    """:param ref: list of tuples as returned by io.load_multi_fasta_gz"""
    l = self.chromosome_lengths = [len(seq[1]) for seq in ref]
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
    parent_list = [p for s in [s1, s2] for p in s.parents]
    for n in range(self.incest_generations):
      if len(set(parent_list)) < len(parent_list):
        return True
      parent_list = [p for s in parent_list for p in s.parents]
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
  def __init__(self, child_factor=2.0):
    """:param child_factor: Average number of children for average fitness parents"""
    self.child_factor = child_factor

  def breed(self, mate_list, generation):
    """
    :param mate_list: list of pairs of parents from which we generate children
    :returns: list of Samples"""
    children = []
    cnt = 0
    for p1, p2 in mate_list:
      avg_fit = (p1.fitness + p2.fitness) / 2.0
      n_children = int(self.child_factor * (avg_fit + 1.0))
      these_children = []
      for n in range(n_children):
        this_child = Sample(generation, cnt)
        this_child.parents = [p1, p2]
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