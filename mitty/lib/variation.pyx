"""This module defines a structure to carry variation information and some functions to perform operations on
genomes (collections of variations). For a description of design choices and algorithms please see the developer
documentation
"""
# Types of het
ABSENT = 0
HOMOZYGOUS = 3
HET_01 = 1
HET_10 = 2
#       00     01     10     11
GT = ['0/0', '0/1', '1/0', '1/1']  # This needs to match rev 3 definitions


cdef class VariationData:
  """Carries POS, stop REF and ALT information."""
  cdef public:
    long int POS, stop
    str REF, ALT

  def __cinit__(self, int _pos=0, int _stop=0, str _ref='', str _alt=''):
    self.POS, self.stop, self.REF, self.ALT = _pos, _stop, _ref, _alt


def new_variation(int _pos=0, int _stop=0, str _ref='', str _alt='', char het=0, char recessive=0, float fitness=0):
  return Variation(VariationData(_pos, _stop, _ref, _alt), het, recessive, fitness)


cdef class Variation:
  """Carries *reference* to a VariationData instance and heterozygosity, recessiveness and fitness data."""
  cdef public:
    VariationData vd
    unsigned char het, recessive, fitness

  # http://docs.cython.org/src/userguide/special_methods.html
  def __cinit__(self, VariationData _vd, unsigned char _het, unsigned char _rec, unsigned char _fit):
    self.vd, self.het, self.recessive, self.fitness = _vd, _het, _rec, _fit

  cpdef bint eq(self, Variation other):
    # We don't consider fitness or recessiveness
    return self.vd.POS == other.vd.POS and self.vd.stop == other.vd.stop and \
           self.vd.REF == other.vd.REF and self.vd.ALT == other.vd.ALT and \
           self.het == other.het

  # http://docs.cython.org/src/userguide/special_methods.html
  def __richcmp__(self, Variation other, int op):
    if op == 2:
      return self.eq(other)
    elif op == 3:
      return not self.eq(other)
    else:
      return False

  def __repr__(self):
    """Return a human friendly printed representation of the variant."""
    return '(POS={0},stop={1},REF={2},ALT={3},het={4},r={5},fit={6})'.\
      format(self.vd.POS, self.vd.stop, self.vd.REF, self.vd.ALT, GT[self.het], self.recessive, self.fitness)


def copy_variant_sequence(c1):
  """c1 - a sequence of variants that need to be copied over. Note that we always share the VariantData"""
  return [Variation(v1.vd, v1.het, v1.recessive, v1.fitness) for v1 in c1]


cdef inline Variation copy_variant(v1):
  return Variation(v1.vd, v1.het, v1.recessive, v1.fitness)


def merge_variants(c1, c2):
  """
  Given an existing chromosome (list of variants) merge a new list of variants into it in zipper fashion

  Args:
    c1 (variant list): The original variants
    c2 (variant list): The proposed new variants

  Returns:
    c3 (variant list): The resultant variant list with collisions arbitrated

  Algorithm:
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
  return c_merge_variants(c1, c2)


cdef inline bint overlap(Variation x, Variation y):
  # Returns true if the footprints of the variations overlap.
  if x is None or y is None: return False
  if y.vd.POS - 1 <= x.vd.POS <= y.vd.stop + 1 or y.vd.POS - 1 <= x.vd.stop <= y.vd.stop + 1 or \
                  x.vd.POS <= y.vd.POS - 1 <= y.vd.stop + 1 <= x.vd.stop:
    # Potential overlap
    if x.het == y.het or x.het == HOMOZYGOUS or y.het == HOMOZYGOUS:
      return True
  return False


#TODO: refactor this to be faster? Not use copy_variant? Code can be made cleaner
cdef c_merge_variants(c1, c2):
  c1_iter, c2_iter = c1.__iter__(), c2.__iter__()
  c3 = []  #deque()
  append = c3.append

  last_new = None
  # Try the zipper
  cdef Variation existing = next(c1_iter, None), denovo = next(c2_iter, None)
  while existing is not None and denovo is not None:
    if overlap(existing, denovo):
      # This will collide, resolve in favor of existing and advance both lists
      append(copy_variant(existing))
      last_new = existing
      existing, denovo = next(c1_iter, None), next(c2_iter, None)
    else:
      if existing.vd.POS <= denovo.vd.POS:  # Zip-in existing
        append(copy_variant(existing))
        last_new = existing
        existing = next(c1_iter, None)
      else:  # Can we zip-in denovo?
        if not overlap(last_new, denovo):
          append(copy_variant(denovo))
          last_new = denovo
        denovo = next(c2_iter, None)  # In either case, we need to advance denovo

  # Now pick up any slack
  if existing is not None:  # Smooth sailing, just copy over the rest
    while existing is not None:
      append(copy_variant(existing))
      existing = next(c1_iter, None)
  else:  # Need to test for overlap before copying over
    while denovo is not None:
      if not overlap(last_new, denovo):
        append(copy_variant(denovo))
        last_new = denovo
      denovo = next(c2_iter, None)  # In either case, we need to advance denovo

  return c3