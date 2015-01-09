#!python
"""This module implements functions for simulating populations of genomes. When run as a script it will produce sets of
VCF files one for each sample in a population and an ancestral_vcf which is the master set of ancestral variations.

Commandline::

  Usage:
    population --pfile=PFILE  [-v]
    population models
    population explain <model_name>

  Options:
    --pfile=PFILE           Name for parameter file
    -v                      Dump detailed logger messages
    models                  List the available denovo models
    explain                 Explain details about the indicated plugin
    <model_name>            The model to explain. If none, explains the parameter file format

Parameter file example::

  {
    "files": {
      "genome": "/Users/kghose/Data/hg38",  # An absolute path is left as is
      "output_dir": "Out",                  # A relative path is taken relative to the location of the *script*
      "output_vcf_prefix": "pop"            # The prefix added to each genome file
                                            # The file names will be pop_sX and the sample names will be sX
                                            # where X 0,1,2 ....
                                            # In addition, there will be a file pop_ancestral_pool.vcf.gz which lists
                                            # The variants used in the ancestral pool
    },
    "rng": {
      "master_seed": 1
    },
    "founder_population": {
      "size": 20,              # Size of population we want
      "p_ancestral": 0.001,    # Probability that a variant in the ancestral pool ends up in a sample
                               # in general p_ancestral * p_model = p_sample
      "p_ancestral_het": 0.9,   # Probability that ancestral variant is het
      "ancestral_models": [        # The list of ancestral variant models should come under this key
        {
          "snp": {                 # name of the model.
            "chromosome": [1, 2],  # Chromosomes to apply this model to
            "phet": 0.5,           # Parameters required by the model. Note that zygosity is ignored
            "p": 0.01
          }
        },
        {                          # We can chain as many models as we wish
          "delete" : {             # We can repeat models if we want
            "chromosome": [1],
            "phet": 0.5,
            "p": 0.01
          }
        }
      ],
      "denovo_models": [           # The list of denovo models should come under this key
        {
          "snp": {                 # name of the model.
            "chromosome": [1, 2],  # Chromosomes to apply this model to
            "phet": 0.5,           # Parameters required by the model
            "p": 0.01
          }
        },
        {                          # We can chain as many models as we wish
          "delete" : {             # We can repeat models if we want
            "chromosome": [1],
            "phet": 0.5,
            "p": 0.01
          }
        }
      ],
    }
  }
"""
__version__ = '1.0.0'

import docopt
import mitty.lib
import mitty.lib.io
import mitty.lib.util
import mitty.lib.genome
from mitty.lib.variation import copy_variant_sequence
import mitty.denovo
import logging
logger = logging.getLogger(__name__)


def founder_population(ancestral_models=[], denovo_models=[], p_a=.95, p_het=0.5, ref=None, pop_size=10, master_seed=1):
  """Generate a founder population

  :param list ancestral_models: List of ancestral models. This will only be run once, so make sure the pool of variants is dense enough for your population size.
  :param list denovo_models: List of denovo models. Will be called for each sample in the population
  :param float p_a: the probability that any given variant from the ancestral pool will end up in a sample
  :param float p_het: the probability that ancestral variants are heterozygous
  :param genome ref: reference genome
  :param int pop_size: size of the founder population
  :param int master_seed: master seed for simulations
  :returns list, genome: list of genomes, ancestral genome (dense list of variants)
  """
  ancestral_model_seed, denovo_model_seed_of_seeds, other_seed = mitty.lib.get_seeds(master_seed=master_seed, size=3)
  ancestral_pool = mitty.denovo.create_variant_list_from_models(ancestral_models, ref, ancestral_model_seed)
  size_of_ancestral_pool = {k: len(g) for k,g in ancestral_pool.iteritems()}

  ancestral_comb_rng, ancestral_het_rng, ancestral_copy_rng = mitty.lib.util.initialize_rngs(master_seed=other_seed, n_rngs=3)
  denovo_model_seeds = mitty.lib.get_seeds(master_seed=denovo_model_seed_of_seeds, size=pop_size)

  founder_pool = []
  for sn in range(pop_size):
    ancestor = {}
    for k, v in ancestral_pool.iteritems():  # Iterate through chromosomes
      ancestral_comb = mitty.lib.util.place_poisson(ancestral_comb_rng, p_a, size_of_ancestral_pool[k])
      ancestral_het = mitty.lib.util.zygosity(num_vars=ancestral_comb.size, phet=p_het, het_rng=ancestral_het_rng, copy_rng=ancestral_copy_rng)
      ancestor[k] = list(copy_variant_sequence(ancestral_pool[k], ancestral_comb, ancestral_het))
    denovo = mitty.denovo.create_variant_list_from_models(denovo_models, ref, denovo_model_seeds[sn])
    founder_pool += [mitty.lib.variation.merge_genomes(ancestor, denovo)]

  return founder_pool, ancestral_pool


def save_population(pop, vcf_prefix='pop'):
  """Given a population of genomes same them as a collection of vcf files with a prefix and sequential numbering."""
  for n, p in enumerate(pop):
    mitty.lib.io.vcf_save_gz(p, vcf_gz_name='{:s}_s{:d}.vcf.gz'.format(vcf_prefix, n), sample_name='s{:d}'.format(n))


def cli():
  """Serves as entry point for scripts"""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    cmd_args = docopt.docopt(__doc__, version=__version__)
  if cmd_args['models']:
    mitty.denovo.print_model_list()
    exit(0)
  if cmd_args['explain']:
    mitty.denovo.explain_model(cmd_args['<model_name>'])
    exit(0)

  level = logging.DEBUG if cmd_args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  import os
  import json
  base_dir = os.path.dirname(cmd_args['--pfile'])     # Other files will be with respect to this
  params = json.load(open(cmd_args['--pfile'], 'r'))

  master_seed = params['rng']['master_seed']

  ref_genome = mitty.lib.genome.FastaGenome(seq_dir=mitty.lib.rpath(base_dir, params['files']['genome']), persist=True)
  ancestral_models = mitty.denovo.load_variant_model_list(params['founder_population']['ancestral_models'])
  denovo_models = mitty.denovo.load_variant_model_list(params['founder_population']['denovo_models'])

  p_a = params['founder_population']['p_ancestral']
  p_het = params['founder_population']['p_ancestral_het']
  pop_size = params['founder_population']['size']
  pop, pool = founder_population(ancestral_models=ancestral_models, denovo_models=denovo_models,
                                 p_a=p_a, p_het=p_het, ref=ref_genome, pop_size=pop_size, master_seed=master_seed)

  output_dir = mitty.lib.rpath(base_dir, params['files']['output_dir'])
  vcf_prefix = params['files']['output_vcf_prefix']
  mitty.lib.io.vcf_save_gz(pool, os.path.join(output_dir, vcf_prefix + '_ancestral_pool.vcf.gz'))
  save_population(pop, vcf_prefix=os.path.join(output_dir, vcf_prefix))


if __name__ == "__main__":
  cli()

#
#
# def chrom_crossover(c1, crossover_idx):
#   """cross_over_idx is a list the same size as c1 with a 1 (indicating a crossover should be done) or 0 (no crossover).
#   """
#   c2 = []
#   for c, idx in zip(c1, crossover_idx):
#     new_c = vcopy(c)  # Valid for no crossover or HOMOZYGOUS
#     if idx == 1 and c.zygosity != HOMOZYGOUS:
#       if c.zygosity == HET_10:
#         new_c.zygosity = HET_01
#       else:
#         new_c.zygosity = HET_10
#     c2 += [new_c]
#   return c2
#
#
# def crossover_event(g1, crossover_idx):
#   return {chrom: chrom_crossover(g1[chrom], crossover_idx[chrom]) for chrom in g1}
#
#
# def pair_one_chrom(c1, c2, which_copy):
#   """Has a resemblance to merge sort BUT WITH GENOMIC VARIANTS!
#   Parameters
#   ----------
#   c1, c2     : tuple list
#                (POS, POS2, REF, ALT, HET)
#   which_copy : tuple
#                e.g. (0, 1) telling us which chromosome of each individual to take
#
#   Returns
#   -------
#   c3         : tuple list
#                (POS, POS2, REF, ALT, HET)
#   """
#   c1_iter, c2_iter = c1.__iter__(), c2.__iter__()
#   c3 = []
#   # Try the zipper
#   l1, l2 = next(c1_iter, None), next(c2_iter, None)
#   while l1 is not None and l2 is not None:
#
#     if (l1.zygosity == HET_10 and which_copy[0] == 1) or (l1.zygosity == HET_01 and which_copy[0] == 0):
#       l1 = next(c1_iter, None)
#       continue
#
#     if (l2.zygosity == HET_10 and which_copy[1] == 1) or (l2.zygosity == HET_01 and which_copy[1] == 0):
#       l2 = next(c2_iter, None)
#       continue
#
#     if vcopy(l1, zygosity=HOMOZYGOUS) == vcopy(l2, zygosity=HOMOZYGOUS):  # Homozygous
#       c3 += [vcopy(l1, zygosity=HOMOZYGOUS)]
#       l1, l2 = next(c1_iter, None), next(c2_iter, None)
#       continue
#
#     if l1.POS <= l2.POS:
#       c3 += [vcopy(l1, zygosity=HET_10)]
#       l1 = next(c1_iter, None)
#     else:
#       c3 += [vcopy(l2, zygosity=HET_01)]
#       l2 = next(c2_iter, None)
#
#   # Now pick up any slack
#   while l1 is not None:
#     if (l1.zygosity == HOMOZYGOUS) or (l1.zygosity == HET_10 and which_copy[0] == 0) or (l1.zygosity == HET_01 and which_copy[0] == 1):
#       c3 += [vcopy(l1, zygosity=HET_10)]
#     l1 = next(c1_iter, None)
#
#   while l2 is not None:
#     if (l2.zygosity == HOMOZYGOUS) or (l2.zygosity == HET_10 and which_copy[1] == 0) or (l2.zygosity == HET_01 and which_copy[1] == 1):
#       c3 += [vcopy(l2, zygosity=HET_01)]
#     l2 = next(c2_iter, None)
#
#   return c3
#
#
# def fertilize_one(g1, g2, which_copy):
#   #return {chrom: pair_one_chrom(c1, c2, wc) for (chrom, c1), (_, c2), wc in zip(g1.iteritems(), g2.iteritems(), which_copy)}
#   return {chrom: pair_one_chrom(g1[chrom], g2[chrom], which_copy[chrom]) for chrom in g1}
#
#
# def place_crossovers_on_chrom(c1, hot_spots, rng):
#   """Stock cross over generator. Place segments based on gaussian distribution around hot_spots
#   The model is a mixture-gaussian model: the hotspot has an amplitude and a width (sd of the gaussian). The amplitude of
#   the gaussian at a certain locus determines the probability of the variants there crossing over. The probability of
#   crossover at any point is the sum of all Gaussians clipped at +1.0
#
#   Parameters
#   ----------
#   c1               : list of tuples
#                      [(start, stop, <any other payload>) ...]. Should be sorted
#   hot_spots        : numpy array
#                      n x 3
#                      [(center, max_p, width) ... ] 0 <= max_p <= 1
#   rng              : random number generator
#
#   Returns
#   -------
#   idx              : numpy array
#                      As used by chrom_crossover
#   """
#   if hot_spots.shape[0] == 0:
#     return [0] * len(c1)  # No hotspots, no crossovers
#   #x = numpy.array(c1, dtype=[('st', float), ('a', 'c'), ('b', 'c'), ('c', 'c'), ('d', 'c')])['st']
#   x = [v.POS for v in c1]
#   X, C = numpy.meshgrid(x, hot_spots[:, 0])
#   _, A = numpy.meshgrid(x, hot_spots[:, 1])
#   _, W = numpy.meshgrid(x, hot_spots[:, 2])
#   p = (A * numpy.exp(-(((X - C) / W) ** 2))).sum(axis=0).clip(0, 1.0)
#   #return numpy.nonzero(rng.uniform(size=p.size) < p)[0]
#   return rng.uniform(size=p.size) < p
#
#
# def place_crossovers(g1, hot_spots, rng):
#   return {chrom: place_crossovers_on_chrom(g1[chrom], hot_spots[chrom], rng) for chrom in g1}
#
#
# def which_copies(g1, rng):
#   """Simple 50/50 coin toss to decide which chromosome copy gets picked during fertilization."""
#   return {chrom: wc for chrom, wc in zip(g1, rng.randint(0, 2, size=(len(g1), 2)))}
#
#
# def get_rngs(seed):
#   """Creates independent RNGs for our random processes from the master seed given."""
#   rng_names = ['cross_over', 'chrom_copy', 'couple_chose']
#   return {k: numpy.random.RandomState(seed=sub_seed)
#           for k,sub_seed in zip(rng_names, numpy.random.RandomState(seed=seed).randint(100000000, size=len(rng_names)))}
#
#
# def spawn(g1, g2, hot_spots={}, rngs={}, num_children=2, ref=None, models=[]):
#   """If you don't want denovo mutations don't send in ref_fp and models"""
#   if g1.keys() != g2.keys():
#     raise RuntimeError('Two genomes have unequal chromosomes')
#   children = []
#   for _ in range(num_children):
#     g1_cross = crossover_event(g1, place_crossovers(g1, hot_spots, rngs['cross_over']))
#     g2_cross = crossover_event(g2, place_crossovers(g2, hot_spots, rngs['cross_over']))
#     if ref is not None and len(models) > 0:  # We want denovo mutations
#       g1_cross = mitty.denovo.apply_variant_models_to_genome(g1_cross, ref, models)
#       g2_cross = mitty.denovo.apply_variant_models_to_genome(g2_cross, ref, models)
#     children.append(fertilize_one(g1_cross, g2_cross, which_copies(g1, rngs['chrom_copy'])))
#   return children
#
#
# # def find_mates(pop, related_level=0, parent_list=[]):
# #   """Related level is the number of generations we need to go back to ensure they are not related."""
# #   our_pop = list(pop)
# #   pairs = []
# #   n = 0
# #   while n < len(our_pop):
# #     if our_pop[n] is None:
# #       n += 1
# #       continue
# #     m = 0
# #     while m < len(pop):
# #
# #     this_pair = [n]
#
#
# def one_generation(pop, hot_spots={}, rngs={}, num_children_per_couple=2, ref=None, models=[]):
#   """We take pairs from the population without replacement, cross over the chromosomes, sprinkle in denovo mutations
#   then shuffle the chromosomes during fertilization. If the population size is odd, we discard one parent"""
#   assert len(pop) >= 2, 'Need at least two parents'
#   mates = rngs['couple_chose'].permutation(len(pop))
#   children = []
#   parents = []
#   for n in range(2*(len(mates)/2))[::2]:  #Todo fix handling of odd sizes
#     children += spawn(pop[mates[n]], pop[mates[n + 1]], hot_spots=hot_spots, rngs=rngs, num_children=num_children_per_couple,
#                       ref=ref, models=models)
#     parents += [(mates[n], mates[n + 1])]
#   return children, parents
#
#
# def population_simulation(ref, denovo_models=[], initial_size=10,
#                           hot_spots={}, rngs={}, num_children_per_couple=2, ss_models=[],
#                           num_generations=10,
#                           store_all_generations=True,
#                           vcf_prefix=''):
#   # The initial population
#   pop = de_novo_population(ref, denovo_models, size=initial_size)
#   logger.debug('Created de novo population of {:d} individuals'.format(len(pop)))
#   save_one_generation(pop, generation=0, vcf_prefix=vcf_prefix)
#
#   parent_list = []
#   # The generational loop
#   this_gen = children = pop
#   for n in range(1, num_generations + 1):
#     children, parents = one_generation(this_gen, hot_spots=hot_spots, rngs=rngs,
#                                        num_children_per_couple=num_children_per_couple,
#                                        ref=ref, models=ss_models)
#     logger.debug('Created generation {:d} ({:d})'.format(n, len(children)))
#     parent_list.append(parents)
#     if store_all_generations:
#       save_one_generation(children, generation=n, vcf_prefix=vcf_prefix)
#     this_gen = children
#
#   if not store_all_generations:
#       save_one_generation(pop, generation=n, vcf_prefix=vcf_prefix)
#
#   return parent_list, children
#
#
# def save_one_generation(pop, generation=0, vcf_prefix='pop'):
#   """Given a population of genomes same them as a collection of vcf files with a prefix and sequential numbering."""
#   for n, p in enumerate(pop):
#     vcf_save_gz(p, vcf_gz_name='{:s}_g{:d}_p{:d}.vcf.gz'.format(vcf_prefix, generation, n),
#                 sample_name='g{:d}p{:d}'.format(generation, n))
#
#
# def save_multiple_generations(gen_pop, vcf_prefix='pop'):
#   for gen, pop in enumerate(gen_pop):
#     save_one_generation(pop, generation=gen, vcf_prefix=vcf_prefix)
#
#
# def main(ref_fname='../../Data/ashbya_gossypii/ashbya.h5', hot_spots={}, params_json={}, master_seed=None,
#          initial_size=10, num_children_per_couple=2, num_generations=20,
#          store_all_generations=True,
#          out_dir='out',
#          out_prefix='pop'):
#   import mitty.denovo as denovo
#   import shutil
#   import os
#   shutil.rmtree(out_dir)
#   os.makedirs(out_dir)
#
#   models = denovo.load_variant_models(params_json['denovo_variant_models'])
#   if master_seed is not None:
#     denovo.apply_master_seed(models, master_seed=int(master_seed))
#     rngs = get_rngs(seed=int(master_seed))
#   else:
#     rngs = get_rngs(seed=1)
#   ss_models = denovo.load_variant_models(params_json['ss_variant_models'])
#   ref = genome.FastaGenome(ref_fname)
#   parent_list, _ = \
#     population_simulation(ref, denovo_models=models, initial_size=initial_size,
#                           hot_spots=hot_spots, rngs=rngs, num_children_per_couple=num_children_per_couple,
#                           ss_models=ss_models,
#                           num_generations=num_generations,
#                           store_all_generations=store_all_generations,
#                           vcf_prefix=os.path.join(out_dir, out_prefix))
#
# if __name__ == "__main__":
#   import json
#   import docopt
#
#   if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
#     docopt.docopt(__doc__, ['-h'])
#   else:
#     args = docopt.docopt(__doc__, version=__version__)
#   if args['plugins']:
#     mitty.denovo.print_plugin_list()
#     exit(0)
#   if args['explain']:
#     mitty.denovo.explain_plugin(args['<plugin>'])
#     exit(0)
#
#   logging.basicConfig(level=logging.WARNING)
#   if args['-v']:
#     logger.setLevel(logging.DEBUG)
#
#   hs = json.load(open(args['--hs'], 'r'))
#   hot_spots = {int(c): numpy.array(v) for c, v in hs.iteritems()}
#   pj = json.load(open(args['--paramfile'], 'r'))
#   main(ref_fname=args['--wg'], hot_spots=hot_spots, params_json=pj, master_seed=args['--master_seed'],
#        initial_size=int(args['--init_size']), num_children_per_couple=int(args['--children']),
#        num_generations=int(args['--gens']), store_all_generations=False,
#        out_dir=args['--outdir'], out_prefix=args['--outprefix'])