"""Example read profile file for reads program

TODO: Potential problem - too much computation offloaded here?

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
paired = True
read_len = 30  # Shorter than this samtools can't align properly
template_len = 250
coverage = 50

#
# def read_params(n=1000, paired=False, rng=None):
#   """Example read_params function for reads.py. This function should return read lengths and template lengths
#   Inputs:
#     n        - how many reads to generate
#     paired   - are we returning paired reads. If not the template len is None
#     rng      - random number generator from reads.py (rng = numpy.random.RandomState(seed))
#   Output:
#     read_len - read lengths. list of tuples. If paired reads then tuples have two elements, otherwise one
#     template_len - length of template. List. None if not paired
#
#   Notes:
#   1. Using the passed random number generator will ensure repeatability wrt the original seed passed to reads.py
#   """
#   mean_rl = 100
#   max_rl = 150
#   rl = rng.poisson(lam=mean_rl, size=n)
#   idx = pylab.find(rl > max_rl)
#   while idx.size > 0:
#     rl2 = rng.poisson(lam=mean_rl, size=idx.size)
#     rl[idx] = rl2
#     idx = pylab.find(rl > max_rl)
#   return rl
#
#
#
#
#
# def read_len(n=1000, rng=None):
#   """Example read length function for reads.py
#   Inputs:
#     n    - how many reads to generate
#     rng  - random number generator from reads.py (rng = numpy.random.RandomState(seed))
#   Output:
#     numpy array of ints.
#     (Any list like element (that can be indexed) is fine)
#
#   Notes:
#   1. Using the passed random number generator will ensure repeatability wrt the original seed passed to reads.py
#   """
#   mean_rl = 100
#   max_rl = 150
#   rl = rng.poisson(lam=mean_rl, size=n)
#   idx = pylab.find(rl > max_rl)
#   while idx.size > 0:
#     rl2 = rng.poisson(lam=mean_rl, size=idx.size)
#     rl[idx] = rl2
#     idx = pylab.find(rl > max_rl)
#   return rl
#
#
# def template_len(n=1000, rng=None):
#   """Example template length function for reads.py
#   Inputs:
#     n    - how many reads to generate
#     rng  - random number generator from reads.py (rng = numpy.random.RandomState(seed))
#   Output:
#     numpy array of ints.
#     (Any list like element (that can be indexed) is fine)
#
#   Notes:
#   1. Using the passed random number generator will ensure repeatability wrt the original seed passed to reads.py
#   """
#   mean_tl = 1000
#   max_tl = 2000
#   tl = rng.poisson(lam=mean_tl, size=n)
#   idx = pylab.find(tl > max_tl)
#   while idx.size > 0:
#     tl2 = rng.poisson(lam=mean_tl, size=idx.size)
#     tl[idx] = tl2
#     idx = pylab.find(tl > max_tl)
#   return tl
#
#
# def read_errors(ideal_read, rng=None):
#   """Example function to generate read errors for reads.py
#   Inputs:
#     ideal_read      - ideal read (stringlike)
#     rng             - random number generator from reads.py (rng = numpy.random.RandomState(seed))
#   Outputs:
#     corrupted_read  - ideal read with read errors inserted
#     phred           - corresponding Phred scores
#   For now, this is just the null model. Will replace this with useful code soon.
#   """
#   return ideal_read, [50] * len(ideal_read)
