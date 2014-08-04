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
import blist
import numpy

HOMOZYGOUS = 0
HET1 = 1
HET2 = 2


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

    chrom += [(start, stop, ref, alt, het)]

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

  correct_chrom = [
      (1, 2, 'C', 'CAA', HET2),
      (3, 6, 'CAG', 'C', HET1),
      (7, 8, 'G', 'T', HET2),
      (9, 12, 'GTT', '', HET1),
      (13, 16, 'GTT', 'TTG', HOMOZYGOUS)
    ]

  chrom = vcf2chrom(vcf.Reader(fsock=io.BytesIO(vcf_str)))
  assert chrom == correct_chrom


def parse_vcf(vcf_rdr, chrom_list):
  """Given a vcf reader load in all the chromosomes."""
  g1 = {}
  for chrom in chrom_list:
    try:
      g1[chrom] = vcf2chrom(vcf_rdr.fetch(chrom, start=0))
    except KeyError:
      g1[chrom] = None

  return g1


def chrom_crossover(c1, cross_over_idx):
  """cross_over_idx is a list the same size as c1 with a 1 (indicating a crossover should be done) or 0 (no crossover).
  """
  c2 = []
  for c, idx in zip(c1, cross_over_idx):
    if idx == 0:
      c2 += [c]
    elif c[4] == HOMOZYGOUS:
      c2 += [c]
    elif c[4] == HET1:
      c2 += [c[:4] + (HET2,)]
    else:
      c2 += [c[:4] + (HET1,)]
  return c2


def chrom_crossover_test():
  c1 = [
      (1, 2, 'C', 'CAA', HET2),
      (3, 6, 'CAG', 'C', HET1),
      (7, 8, 'G', 'T', HET2),
      (9, 12, 'GTT', '', HET1),
      (13, 16, 'GTT', 'TTG', HOMOZYGOUS),
      (23, 26, 'GTT', 'TTG', HOMOZYGOUS)
    ]
  cross_over_idx = [1, 1, 1, 0, 1, 0]
  correct_chrom = [
      (1, 2, 'C', 'CAA', HET1),
      (3, 6, 'CAG', 'C', HET2),
      (7, 8, 'G', 'T', HET1),
      (9, 12, 'GTT', '', HET1),
      (13, 16, 'GTT', 'TTG', HOMOZYGOUS),
      (23, 26, 'GTT', 'TTG', HOMOZYGOUS)
    ]
  c1_c = chrom_crossover(c1, cross_over_idx)
  assert correct_chrom == c1_c, c1_c
  assert c1[0] == (1, 2, 'C', 'CAA', HET2)  # We shouldn't be changing our original


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

    if (l1[4] == HET1 and which_copy[0] == 1) or (l1[4] == HET2 and which_copy[0] == 0):
      l1 = next(c1_iter, None)
      continue

    if (l2[4] == HET1 and which_copy[1] == 1) or (l2[4] == HET2 and which_copy[1] == 0):
      l2 = next(c2_iter, None)
      continue

    if l1[:4] == l2[:4]:  # Homozygous
      c3 += [l1[:4] + (HOMOZYGOUS,)]
      l1, l2 = next(c1_iter, None), next(c2_iter, None)
      continue

    if l1[0] <= l2[0]:
      c3 += [l1[:4] + (HET1,)]
      l1 = next(c1_iter, None)
    else:
      c3 += [l2[:4] + (HET2,)]
      l2 = next(c2_iter, None)

  # Now pick up any slack
  while l1 is not None:
    c3 += [l1[:4] + (HET1,)]
    l1 = next(c1_iter, None)

  while l2 is not None:
    c3 += [l2[:4] + (HET2,)]
    l2 = next(c2_iter, None)

  return c3


def pair_one_chrom_test():
  c1 = [
    (1, 2, 'C', 'CAA', HET2),  # Tests homozygosity
    (3, 6, 'CAG', 'C', HET1),  # Tests both variants are not on the copies chosen
    (17, 18, 'G', 'T', HET2),  # Test zipper (several variants should come from c2 before we get to this)
    (23, 26, 'GTT', 'TTG', HOMOZYGOUS)  # Tests handling of homozygous variants
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
  ]
  assert correct_pairing == pair_one_chrom(c1, c2, which_copy)


def fertilize_one(g1, g2, which_copy):
  return {chrom: [g1[chrom][p[0]], g2[chrom][p[1]]] for chrom, p in zip(sorted(g1.keys()), which_copy)}


def fertilize_one_test():
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


def place_crossovers(g1, hot_spots, rng):
  """Stock cross over generator. Place segments based on gaussian distribution around hotspots

  Parameters
  ----------
  g1               : list of tuples
                     [(start, stop, <any other payload>) ...]. Should be sorted
  hot_spots        : list
  rng              : random number generator

  Returns
  -------
  idx              : list
                     As used by chrom_crossover
  """
  idx = sorted([i + r for r, i in zip(rng.randn(len(hot_spots)), hot_spots)])
  for i in idx:
    pass


  min_coord = min(cpy1[0][0], cpy2[0][0])
  max_coord = max(cpy1[-1][1], cpy2[-1][1])
  x_over_pts = rng.uniform(low=min_coord, high=max_coord, size=2 * cnt)
  x_over_pts.sort()

  # st, nd = 0
  # while nd <
  #
  # rng.poisson(lam, 2)




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




def pairing(g1, g2, rng, num_children=2):
  if g1.keys() != g2.keys():
    raise RuntimeError('Two genomes have unequal chromosomes')
  return [one_pairing(g1, g2, rng.randint(2, size=(len(g1), 2))) for _ in range(num_children)]






def gen2vcf(g1):
  pass

