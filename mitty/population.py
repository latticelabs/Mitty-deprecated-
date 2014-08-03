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
vcf2chrom(vcf_reader)  -> g2 read in one chromosome from a VCF file
write_to_vcf(fname, g1) write the data to a compressed, indexed VCF file




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
import blist
import numpy


def vcf2chrom(vcf_rdr):
  """Given a vcf reader corresponding to one chromosome, read in the variant descriptions into our format. The result is
  sorted if the vcf file is sorted.
  """
  chrom = [[], []]

  for variant in vcf_rdr:
    alt = variant.ALT[0].sequence if variant.ALT[0] is not None else ''
    ref = variant.REF or ''
    start = variant.POS  # Note, we are in VCF coordinates!
    stop = variant.POS + len(ref)
    v = [(start, stop, ref, alt)]

    try:
      if variant.samples[0].gt_nums[0] == '1':
        chrom[0] += v
    except IndexError:
        chrom[0] += v
    try:
      if variant.samples[0].gt_nums[2] == '1':
        chrom[1] += v
    except IndexError:
        chrom[1] += v

  return chrom


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

  cpy1 = [
    (3, 6, 'CAG', 'C'),
    (9, 12, 'GTT', ''),
    (13, 16, 'GTT', 'TTG')
  ]
  cpy2 = [
    (1, 2, 'C', 'CAA'),
    (7, 8, 'G', 'T'),
    (13, 16, 'GTT', 'TTG')
  ]

  chrom = vcf2chrom(vcf.Reader(fsock=io.BytesIO(vcf_str)))
  assert chrom[0] == cpy1
  assert chrom[1] == cpy2


def chrom_merge(cpy1, cpy2):
  """Same as merge sort BUT ON GENOMIC VARIANTS!
  This returns a sorted list ready to go into a vcf file

  Returns
  -------
  (POS, POS2, REF, ALT, HET)
  """
  c1, c2 = sorted(cpy1), sorted(cpy2)
  len_c1, len_c2 = len(c1), len(c2)
  vcf_lines = []

  p1, p2 = 0, 0

  while p1 < len_c1 and p2 < len_c2:
    if c1[p1] == c2[p2]:  # Homozygous
      vcf_lines.append(c1[p1] + ('1/1',))
      p1 += 1
      p2 += 1
      continue

    if c1[p1][0] <= c2[p2][0]:
      vcf_lines.append(c1[p1] + ('1/0',))
      p1 += 1
    elif c2[p2][0] < c1[p1][0]:
      vcf_lines.append(c2[p2] + ('0/1',))
      p2 += 1

  while p1 < len_c1:
    vcf_lines.append(c1[p1] + ('1/0',))
    p1 += 1

  while p2 < len_c2:
    vcf_lines.append(c2[p2] + ('0/1',))
    p2 += 1

  return vcf_lines


def chrom_merge_test():
  cpy1 = [
    (1, 2, 'C', 'CA'),
    (3, 6, 'CAG', 'C'),
    (13, 16, 'GTT', 'TTG'),
    (9, 12, 'GTT', ''),
    (20, 23, 'GTT', 'TTG'),
  ]
  cpy2 = [
    (20, 21, 'T', 'G'),
    (1, 2, 'C', 'CAA'),
    (13, 16, 'GTT', 'TTG'),
    (7, 8, 'G', 'T')
  ]

  correct_vcf = [
    (1, 2, 'C', 'CA', '1/0'),
    (1, 2, 'C', 'CAA', '0/1'),
    (3, 6, 'CAG', 'C', '1/0'),
    (7, 8, 'G', 'T', '0/1'),
    (9, 12, 'GTT', '', '1/0'),
    (13, 16, 'GTT', 'TTG', '1/1'),
    (20, 23, 'GTT', 'TTG', '1/0'),
    (20, 21, 'T', 'G', '0/1')
  ]

  assert correct_vcf == merge(cpy1, cpy2)


def get_rngs(seed):
  return [numpy.random.RandomState(seed=sub_seed) for sub_seed in numpy.random.RandomState(seed=seed).randint(100000000, size=3)]


def one_pairing(g1, g2, which_copy):
  return {chrom: [g1[chrom][p[0]], g2[chrom][p[1]]] for chrom, p in zip(sorted(g1.keys()), which_copy)}


def one_pairing_test():
  g1 = {
    '1': [[(1, 2, 'C', 'CA')], [(4, 5, 'G', 'T')]],
    '2': [[(3, 6, 'CAG', 'C')], [(8, 9, 'A', 'T')]]
  }
  g2 = {
    '1': [[(2, 5, 'CTT', 'C')], [(8, 11, 'GTT', 'G')]],
    '2': [[(3, 4, 'C', 'T')], [(7, 10, 'GAA', 'G')]]
  }
  g3 = one_pairing(g1, g2, [(0, 0), (1, 1)])
  correct_g3 = {
    '1': [[(1, 2, 'C', 'CA')], [(2, 5, 'CTT', 'C')]],
    '2': [[(8, 9, 'A', 'T')], [(7, 10, 'GAA', 'G')]]
  }
  assert g3 == correct_g3

  g3 = one_pairing(g1, g2, [(0, 1), (1, 0)])
  correct_g3 = {
    '1': [[(1, 2, 'C', 'CA')], [(8, 11, 'GTT', 'G')]],
    '2': [[(8, 9, 'A', 'T')], [(3, 4, 'C', 'T')]]
  }
  assert g3 == correct_g3


def pairing(g1, g2, rng, num_children=2):
  if g1.keys() != g2.keys():
    raise RuntimeError('Two genomes have unequal chromosomes')
  return [one_pairing(g1, g2, rng.randint(2, size=(len(g1), 2))) for _ in range(num_children)]


def chrom_crossover(cpy1, cpy2, cross_over_idx):
  """cross_over_idx is a list of two lists. Each list is the same length as the corresponding chromosome and indicates
  which copy the variant goes to. For example, if there is no cross-over, the list for the first copy will be [0,0,...]
  and the list for the second copy will be [1,1,...]. If the first 2 variants of cpy1 should go to cpy2 but cpy2 is
  otherwise altered we will have [1,1,0,0...] and [1,1,...]

  Returns sorted variants

  """
  all_cpy = cpy1 + cpy2
  all_idx = cross_over_idx[0] + cross_over_idx[1]
  return sorted([va for va, cidx in zip(all_cpy, all_idx) if cidx == 0]), \
         sorted([va for va, cidx in zip(all_cpy, all_idx) if cidx == 1])


def chrom_crossover_test():
  cpy1 = [
    (1, 2, 'C', 'CA'),
    (3, 6, 'CAG', 'C'),
    (9, 12, 'GTT', ''),
    (13, 16, 'GTT', 'TTG'),
    (20, 23, 'GTT', 'TTG')
  ]
  cpy2 = [
    (1, 2, 'C', 'CAA'),
    (7, 8, 'G', 'T'),
    (13, 16, 'GTT', 'TTG'),
    (20, 21, 'T', 'G'),
  ]
  cross_over_idx = [
    [0, 1, 0, 1, 1],
    [0, 1, 0, 1]
  ]

  correct_cpy1 = sorted([
    (1, 2, 'C', 'CA'),
    (9, 12, 'GTT', ''),
    (1, 2, 'C', 'CAA'),
    (13, 16, 'GTT', 'TTG'),
  ])
  correct_cpy2 = sorted([
    (3, 6, 'CAG', 'C'),
    (13, 16, 'GTT', 'TTG'),
    (20, 23, 'GTT', 'TTG'),
    (7, 8, 'G', 'T'),
    (20, 21, 'T', 'G')
  ])
  assert correct_cpy1, correct_cpy2 == chrom_crossover(cpy1, cpy2, cross_over_idx)


def place_crossovers(cpy1, cpy2, rng, cnt=10):
  """Stock cross over generator. Place segments uniformly.

  Parameters
  ----------
  cpy1             : list of tuples
                     [(start, stop, <any other payload>) ...]. Should be sorted
  cpy2             : list of tuples
                     [(start, stop, <any other payload>) ...]. Should be sorted
  rng              : random number generator
  cnt              : int
                     How many cross over segments

  Returns
  -------
  idx              : pair of lists
                     As used by chrom_crossover
  """
  min_coord = min(cpy1[0][0], cpy2[0][0])
  max_coord = max(cpy1[-1][1], cpy2[-1][1])
  x_over_pts = rng.uniform(low=min_coord, high=max_coord, size=2 * cnt)
  x_over_pts.sort()

  st, nd = 0
  while nd <

  rng.poisson(lam, 2)



def gen2vcf(g1):
  pass

