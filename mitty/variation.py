"""This module defines a single, simple, named_tuple that is used in several places by Mitty and should be used by
all the mutation plugins

(start, stop, REF, ALT, het)

"""
from collections import namedtuple
# Types of het
HOMOZYGOUS = 0
HET1 = 1
HET2 = 2
Variation = namedtuple('Variation', 'start, stop, REF, ALT, het')