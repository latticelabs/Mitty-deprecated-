"""This module defines a structure to carry variation information and some functions to perform operations on
genomes (collections of variations). For a description of design choices and algorithms please see the developer
documentation
"""
import itertools  # Because zip does not play well with iterators

# Types of zygosity
cpdef enum:
  ABSENT = 0
  HET_01 = 1
  HET_10 = 2
  HOM = 3
#       00     01     10     11
GT = ['0|0', '0|1', '1|0', '1|1']  # This needs to match rev 3 definitions

# Copy 0 -> 1|0
# Copy 1 -> 0|1

cdef class Variant:
  """Variant(pos, stop, REF, ALT)
  A lightweight class to carry a single data

  Attributes:
    pos    - position of data
    stop   - where does footprint of data on reference (position of last base + 1 in the REF entry)
    hash   - used to check equality between variants
    ref    - reference sequence
    alt    - alternate sequence

  Each Variant is placed in a master list dictionary. The index in the dictionary is computed from the pos and the
  number allele this is at that position. This index is then used as the rowid when we save the data to the database.
  """
  cdef public:
    unsigned long pos, stop, hash, index
    bytes ref, alt

  # def __cinit__(self, unsigned long pos, unsigned long stop, bytes ref, bytes alt):
  #   self.pos, self.stop, self.ref, self.alt = pos, stop, ref, alt
  #   self.hash = self.pos ^ hash(self.ref) ^ hash(self.alt)

  def __repr__(self):
    return '({:d}, {:d}, {:s}, {:s})'.format(self.pos, self.stop, self.ref, self.alt)

  def __richcmp__(self, Variant other, int op):
    if op == 2:
      return other.hash == self.hash
    elif op == 3:
      return not (other.hash == self.hash)

  def as_tuple(self):
    return self.pos, self.stop, self.ref, self.alt


cpdef Variant new_variant(unsigned long pos, unsigned long stop, bytes ref, bytes alt):
  cdef Variant v = Variant()
  v.pos, v.stop, v.ref, v.alt = pos, stop, ref, alt
  v.hash = v.pos ^ hash(v.ref) ^ hash(v.alt)
  return v


cpdef Variant add_novel_variant_to_master(Variant v, dict ml):
  """add_novel_variant_to_master(Variant v, dict ml)
  Add this new variant to the dictionary (master list). Intelligently handle case where we have an existing allele
  at the same locus.

  :param v: This variant in the form of a Variation instance
  :param ml: Python dictionary storing all variants in a chromosome
  :returns: v or existing variant in master list to which this variant is identical
  """
  cdef unsigned short n = 0
  cdef unsigned long index = v.pos << 16 | n
  while index in ml:  # Find us an empty spot
    if ml[index].hash == v.hash:
      return ml[index] # This is identical to an existing variant
    n += 1
    index = v.pos << 16 | n
  ml[index] = v
  v.index = index
  return v


cdef class GTVariant:
  """GTVariant(gt, Variant)

  Represents a data in a sample. It points to the Variant information and carries genotype info. It also has a
  pointer to the next GTVariant so we can chain it together into a Chromosome

  Attributes:
    gt      - zygosity (genotype) information
    data    - (pointer to) the Variant in the master list
  """
  cdef public:
    unsigned char gt
    Variant data
    GTVariant next

  # def __cinit__(self, unsigned char gt, Variant variant):
  #   self.gt, self.data, self.next = gt, variant, None

  def __repr__(self):
    return '{:s} {:s}'.format(self.data, GT[self.gt])

  def __richcmp__(self, GTVariant other, int op):
    if op == 2:
      return other.data.hash == self.data.hash and self.gt == other.gt
    elif op == 3:
      return not (other.data.hash == self.data.hash and self.gt == other.gt)


cpdef GTVariant new_gt_variant(unsigned char gt, Variant variant):
  cdef GTVariant sv = GTVariant()
  sv.gt, sv.data, sv.next = gt, variant, None
  return sv


def create_gtv_iterable(pos, ref, alt, gt):
  """Given lists/generators of the bare variant data (e.g. from a plugin) return us an iterator that will produce a
  stream of individual SampleVariants."""
  for p, r, a, g in itertools.izip(pos, ref, alt, gt):
    yield new_gt_variant(g, new_variant(p, p + len(r), r, a))


cdef class Chromosome:
  """Represents a chromosome as a linked list of GTVariants. We always go sequentially through a chromosome, and for
  denovo variants we may add GTVariants into the middle of the chromosome which is a good fit for a ll

  Rules:
    There is a dummy node, which head (and, initially, tail) point to.
    append: attaches to tail.next and repositions tail. This is how we create the list
    insert: attaches in between cursor and cursor.next
    advance: advances cursor and returns cursor.next

  Note:
    The implementation is specific to our use case, most notably in the use of the 'cursor' pointer and the 'advance'
    method - which returns cursor.next, rather than cursor. Mutable state (in the form of the 'cursor' pointer)
    requires us to remember to rewind when needed. Since this is not a general implementation of a data structure,
    these are acceptable choices.
  """
  cdef public:
    GTVariant head, cursor, tail
    unsigned long length

  def __cinit__(self):
    self.head = self.tail = new_gt_variant(ABSENT, new_variant(0, 0, b'', b''))  # Dummy node
    self.cursor, self.length = None, 0

  def __len__(self):
    return self.length

  def __iter__(self):
    return ChromosomeIterator(self)

  def to_list(self):
    return [sv for sv in self]

  cpdef rewind_cursor(self):
    self.cursor = None

  cpdef GTVariant advance(self):
    if self.cursor is None:  # We are at the head of the list
      self.cursor = self.head
    else:
      if self.cursor.next is not None:  # OK to advance
        self.cursor = self.cursor.next
    return self.cursor.next  # Will be None for empty list or end of list

  cpdef GTVariant last(self):
    if self.head == self.tail:  # Empty list
      return None
    else:
      return self.tail

  cpdef insert(self, GTVariant sv):
    """Place sv right after cursor. cursor will point to original cursor.next after this. In practice, this places
    sv right before the variant advance just spit out"""
    sv.next = self.cursor.next
    self.cursor.next = sv
    self.cursor = sv
    if sv.next is None:
      self.tail = sv
    self.length += 1

  cpdef append(self, GTVariant sv):
    """This is how we grow the list."""
    self.tail.next = sv
    self.tail = sv
    self.length += 1


cdef class ChromosomeIterator:
  """A class that lets us iterate over the sample."""
  cdef GTVariant this

  def __cinit__(self, Chromosome s):
    self.this = s.head.next  # Head is a dummy!

  def __next__(self):
    if self.this is None:
      raise StopIteration()
    else:
      result = self.this
      self.this = self.this.next
      return result


cdef inline bint overlap(GTVariant v1, GTVariant v2):
  """Returns true if the footprints of the variations overlap. This is used when applying denovo mutations to a genome

  :param v1: GTVariant from original sample
  :param v2: Proposed GTVariant
  :returns True/False: Bool
  """
  if v2.data.pos - 1 <= v1.data.pos <= v2.data.stop + 1 or v2.data.pos - 1 <= v1.data.stop <= v2.data.stop + 1 or \
     v1.data.pos <= v2.data.pos - 1 <= v2.data.stop + 1 <= v1.data.stop:  # Potential overlap
    if v1.gt == v2.gt or v1.gt == HOM or v2.gt == HOM:  # Definite overlap
      return True
  return False


cpdef add_denovo_variants_to_sample(Chromosome s, dnv, dict ml):
  """add_denovo_variants_to_sample(s, dnv, ml)
  Given an existing Chromosome, s,  merge a new list of SampleVariants (dnv) into it in zipper fashion. s has
  priority (collisions are resolved in favor of s). As new variants are accepted, add them to the master list.

  :param s: The original variant list
  :param dnv: The proposed new variants (iterator, convenient to use create_gtv_iterable)
  :param ml: master list of Variants

  s and ml are modified in place

  Algorithm::

    o(x,y) = True if x and y overlap
           = False otherwise
    c1 = existing list of variants
    c2 = denovo list of variants
    c3 = new list of variants being built

    o(c1, c2) = True: add(c1) c1++, c2++
              = False:
                c1 < c2 ? add(c1) c1++
                else:
                  o(c3, c2) = True: c2++
                            = False: add(c2), c2++

  """
  s.rewind_cursor()  # One of the evils of mutable state. Is there a way to avoid? Perhaps by making shallow copies?
  cdef:
    GTVariant s1 = s.advance(), s2 = next(dnv, None)
  while s1 is not None and s2 is not None:
    if overlap(s1, s2): # This will collide. Advance dnv and redo. Fixes case for merge_test8
      s2 = next(dnv, None)
    else:  # No collision
      if s1.data.pos <= s2.data.pos:  # Advance s until we come to s2
        s1 = s.advance()
      else:  #s1 is past s2 and there is no collision, good to add
        s2.data = add_novel_variant_to_master(s2.data, ml)  # If this denovo already exists in the master list, use that
        s.insert(s2)
        s2 = next(dnv, None)

  while s2 is not None:  # All the remaining denovo variants can just be appended to the sample
    if s.last() is None or not overlap(s.last(), s2):
      s2.data = add_novel_variant_to_master(s2.data, ml)
      s.append(s2)
    s2 = next(dnv, None)


cdef append_sv_copy(Chromosome s, int child_cp, GTVariant sv, int parent_cp):
  """Append correct copy of GTVariant to s."""
  cdef:
    unsigned char gt
    GTVariant new_sv
  if (sv.gt == HOM) or (sv.gt == HET_01 and parent_cp == 1) or (sv.gt == HET_10 and parent_cp == 0):
    s.append(new_gt_variant(HET_10 if child_cp == 0 else HET_01, sv.data))


cdef append_sv_copies(Chromosome s, GTVariant sv1, int cp1, GTVariant sv2, int cp2):
  """Use this function when both variants have the same pos and we need to arbitrate genotype"""
  # If the variant is absent from one of the copies, don't have to worry about this
  if (sv1.gt == HET_10 and cp1 == 1) or (sv1.gt == HET_01 and cp1 == 0):
    append_sv_copy(s, 1, sv2, cp2)
  elif (sv2.gt == HET_10 and cp2 == 1) or (sv2.gt == HET_01 and cp2 == 0):
    append_sv_copy(s, 0, sv1, cp1)
  else:  # Both copies of the child genome will have a copy of sv1 abd sv2. Now we have to check more deeply
    if sv1.data == sv2.data:  # This will be HOM
      s.append(new_gt_variant(HOM, sv1.data))
    else:  # OK, two variants at the same spot, but they are different
      s.append(new_gt_variant(HET_10, sv1.data))  # They will appear in a certain order
      s.append(new_gt_variant(HET_01, sv2.data))


cdef class CrossOverPoints:
  """We needed this rather than a simple iterator because we have to handle the paradoxes that arise when the cross over
  point is in the middle of a footprint and adjust the cross over accordingly"""
  cdef:
    list point_list
    int current_point, current_copy, current_point_idx, max_point_idx, max_pos

  def __cinit__(self, list point_list, int max_pos, int copy):
    self.point_list, self.max_pos, self.current_copy, self.current_point_idx = point_list, max_pos, copy, 0
    self.max_point_idx = len(self.point_list)
    if self.current_point_idx < self.max_point_idx:
      self.current_point = self.point_list[self.current_point_idx]
    else:
      self.current_point = max_pos

  cdef int update_crossover(self, Variant v):
    """Given the current variant, update the current_copy variable as needed and return that"""
    if v.pos < self.current_point < v.stop:
      # If the cross-over point falls in the middle of a deletion it will lead to a paradox, so we:
      # a) move the cross over point virtually
      self.current_point = v.stop
      # b) Advance the index so that the next advance will point it to beyond this cross over
      while self.current_point_idx < self.max_point_idx:
        if self.point_list[self.current_point_idx] <= self.current_point:
          self.current_point_idx += 1
        else:
          self.current_point_idx -= 1  # We will advance this, so we need to rewind it a bit
          break
      # At the end of this loop, either our next point is past the moved point or we have run out of points
      # No need to flip copies in this case
    else:  # We need to figure out if we have gone past the last cross over and need to advance
      while v.pos > self.current_point:
        self.current_copy = 1 - self.current_copy
        self.current_point_idx += 1
        if self.current_point_idx >= self.max_point_idx:
          self.current_point_idx = self.max_point_idx
          self.current_point = self.max_pos
          break
        else:
          self.current_point = self.point_list[self.current_point_idx]
    return self.current_copy


cpdef Chromosome pair_chromosomes(Chromosome s1, list cross_over1, int chrom_copy1, Chromosome s2, list cross_over2, int chrom_copy2):
  """pair_chromosomes(s1, cross_over_points1, int cp1, s2, cross_over_points2, int cp2)
  Starting with a pair of parent samples, apply crossover points and pair the resulting gametes to make a child genome

  :param s1: first sample genome
  :param cross_over1: list of cross over points
  :param chrom_copy1: which copy of chrom to use for gamete [0, 1]
  :param s2: second sample genome
  :param cross_over2: list of cross over points
  :param chrom_copy2: which copy of chrom to use for gamete [0, 1]
  :returns s3: child genome

  Algorithm:
    * Zip in sample variants from both parents in merge sort fashion until done.
    * For each parent, pick variants from one copy, until we reach a cross over point, then flip the copy
    * By convention, parent0 supplies chrom copy 0 and parent1 supplies chrom copy 1 of the child
  """
  s1.rewind_cursor()  # One of the evils of mutable state. Is there a way to avoid? Perhaps by making shallow copies?
  s2.rewind_cursor()

  cdef:
    int cp1 = chrom_copy1, cp2 = chrom_copy2
    unsigned long co_pt1, co_pt2, p1, p2, st1, st2, max_pos
    Chromosome s3 = Chromosome()
    GTVariant sv1 = s1.advance(), sv2 = s2.advance()

  if s1.last() is not None:
    max_pos = s1.last().data.pos + 1
  else:
    max_pos = 1
  if s2.last() is not None:
    max_pos = max(s2.last().data.pos + 1, max_pos)

  cdef CrossOverPoints cop1 = CrossOverPoints(cross_over1, max_pos, chrom_copy1)
  cdef CrossOverPoints cop2 = CrossOverPoints(cross_over2, max_pos, chrom_copy2)

  while sv1 is not None or sv2 is not None:
    p1 = sv1.data.pos if sv1 else max_pos
    p2 = sv2.data.pos if sv2 else max_pos
    cp1 = cop1.update_crossover(sv1.data)
    cp2 = cop2.update_crossover(sv2.data)

    if p1 < p2:    # Nothing tricky, just add sv1 to the sample
      append_sv_copy(s3, 0, sv1, cp1)  # By convention, parent1 contributes to copy0 of child chromosome
      sv1 = s1.advance()
    elif p2 < p1:  # Nothing tricky, just add sv2 to the sample
      append_sv_copy(s3, 1, sv2, cp2)  # By convention, parent2 contributes to copy1 of child chromosome
      sv2 = s2.advance()
    else:  # A little bit tricky - need to check to see if we end up with a HOM variant, or two HETs
      append_sv_copies(s3, sv1, cp1, sv2, cp2)
      sv1 = s1.advance()
      sv2 = s2.advance()

  return s3





# cpdef merge_variants(list c1, list c2):
#   """merge_variants(c1, c2)
#   Given an existing chromosome (list of variants) merge a new list of variants into it in zipper fashion. c1 has
#   priority (collisions are resolved in favor of c1)
#
#   :param c1: The original variants (iterable)
#   :param c2: The proposed new variants (iterable)
#   :rtype list: The resultant data list with collisions arbitrated
#
#   Algorithm::
#
#     o(x,y) = True if x and y overlap
#            = False otherwise
#     c1 = existing list of variants
#     c2 = denovo list of variants
#     c3 = new list of variants being built
#
#     o(c1, c2) = True: add(c1) c1++, c2++
#               = False:
#                 c1 < c2 ? add(c1) c1++
#                 else:
#                   o(c3, c2) = True: c2++
#                             = False: add(c2), c2++
#
#   """
#   cdef:
#     int n1_max = len(c1), n2_max = len(c2)
#     int n1 = 0, n2 = 0, n3 = 0
#     Variant v1, v2
#     Genotype g1, g2
#     list c3 = [None] * (n1_max + n2_max)  # This is the maximum size of c3
#   while n1 < n1_max and n2 < n2_max:
#     g1, g2 = c1[n1], c2[n2]
#     v1, v2 = ml[g1.index], ml[g2.index]
#     if overlap(v1, v2): # This will collide. Advance c2 and redo. Fixes case for merge_test8
#       n2 += 1
#     else:  # No collision
#       if v1.pos <= v2.pos:  # Zip-in c1 as it comes before c2
#         c3[n3] = g1
#         n1 += 1
#         n3 += 1
#       else:  # Zip-in c2 as it comes before next c1
#         if n3==0 or not overlap(c3[n3 - 1], v2):
#           c3[n3] = v2
#           n3 += 1
#         n2 +=1  # Need to advance c2 anyway
#
#   # Now copy over slack
#   while n1 < n1_max:  # Smooth sailing, just copy over the rest of the original (c1)
#     c3[n3] = v1
#     n1 += 1
#     n3 += 1
#
#   while n2 < n2_max:  # Need to test each new (c2) for clashes with itself
#     if n3==0 or not overlap(c3[n3 - 1], v2):
#       c3[n3] = c2[n2]
#       n3 += 1
#     n2 += 1  # Need to advance c2 anyway
#
#   return c3[:n3]


# def copy_missing_chromosomes(dict g1, dict g2):
#   """Copy any chromosomes found in g2 but not in g1 onto g1. g1 is changed in place
#
#   :param dict g1: genome
#   :param dict g2: genome
#   :returns: changes g1 in place"""
#   cdef list missing_chrom = set(g2.keys()) - set(g1.keys())
#   cdef int ch
#   for ch in missing_chrom:
#     g1[ch] = list(g2[ch])
#
#
# def merge_genomes(dict g1, dict g2):
#   """Given two genomes run merge_variants(c1, c2) on each chromosome. g1 has priority. In the end copy any chromosomes
#   present in g2 but not in g1 into g1 (so that g1 is a superset of g2)
#
#   :param dict g1: genome
#   :param dict g2: genome
#   :returns: new genome
#   """
#   cdef dict g3 = {}
#   for chrom, c2 in g2.iteritems():
#     g3[chrom] = merge_variants(g1.get(chrom, []), c2)
#   copy_missing_chromosomes(g3, g1)  # The previous loop merges all chr in c2.
#                                     # Now, we need to consider all chr in c1 but not in c2
#   return g3