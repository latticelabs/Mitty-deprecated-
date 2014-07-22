The mechanics of simulating variations
=====================================

Scheme
------
* Each variant is described as a sequence of substitutions, insertions or deletions
* substitutions don't change indexing, deletions remove indexes, insertions cause repeat indexes (POS)
* Each variant has one footprint
* Use a mask to arbitrate collisions.
* 1 base buffer between variants
* list of variants in arbitrary order

Variant description:

higher order description  ()  - as needed to describe complex variants - will make VCF more sophisticated as needed later)
footprints                (het, chrom, pos_st, pos_nd) ... used for collision detection
vcf list                  (chrom, pos, id, ref, alt, qual, filter, info, format, sample) ... (as many as needed) - for simple VCF file

het  1 -> variant on copy 1
     2 -> variant on copy 2
     3 -> variant on both copies

Each plugin mass produces variants on a given chromosome and returns lists of tuples as described above. 

The footprint is used to arbitrate mutations and discard those that collide. 

Those that pass are written to the VCF file and to the operations pkl file.

vcf2seq
When we are mutating the sequence, our rules are:

if ref = '.' it's an insertion at pos                        index is repeated
if alt = '.' it's a deletion starting at pos                 index is deleted
if length(ref) = length(alt) = 1 it's a SNP (substitution)          indexes are unchanged
if length(ref) = length(alt) > 1 it's an inversion (substitution)
 
Reconstruction
--------------
Assume the VCF is sorted (which tabix will do for us) 











foot_chrom  - foot print (check for collisions)
foot_start
foot_stop
description - english description of variant
vcf         - ( ) ... 
operations  - (op_name, seq, len, chrom, pos) ....
      
      
   _ type name
  /
(type, het, footprint, )
        \            \
         \            \_ list of variant description tuples
          \
           \_ 





Design choices and explanation
------------------------------
There are a wide variety of mutations that happen in the human genome and they range in size from changes in single bases
to deletions or insertions of 200kb lengths sequences.

To computationally describe these mutations we only require three operations 

  - single base modification
  - insertion
  - deletion
  
The single base modification is not strictly needed - it is a combination of deletion followed by insertion - but it
makes some things simpler for us.

Using combinations these operations allows us to describe all biologically defined variants. For example
  
    SNP                -  substitution
    Insert             -  insertion
    Deletion           -  deletion
    Inversion          -  substitution
    Tandem repeat      -  insertion (s)
    Mobile element     -  deletion + insertion


When we simulate variations there are two things we wish to keep in mind

1. We don't allow variations to overlap
1. There should be at least one [buffer base](#bufferbases) between variants

If either of these rules were to be violated the resulting variations in the sequence could have non unique descriptions
which make it difficult for our purposes of algorithm testing.
    
Describing variants and avoiding collisions
-------------------------------------------
When we generate variants we want to have an unambiguous description of the variation both in terms of the VCF (which
the rest of the world uses) and in terms of our three elementary sequence operations, which we will be using internally. 

Variants may need to be described by operations on different chromosomes, so our coordinate system needs to include the 
chromosome.
 
When we generate a variant, we need to make sure it does not violate our two collision rules above. To ensure this we
use a position mask and the notion of a footprint. A footprint is a section of the reference sequence that is blocked
off by an existing variant. We can not place any other variation here as this would invalidate existing variation 
descriptions.

Some common footprints are described here:

    SNP                -  single base modification
    Insert             -  insertion
    Deletion           -  deletion
    Inversion          -  deletion, followed by insertion
    Tandem repeat      -  insertion (s)
    Mobile element     -  deletion + insertion



Generating the sequence
-----------------------
A simple way of generating the variant sequence is to arrange all the variants by position, start from the beginning of
the sequence and copy over invariant parts. When a variation is encountered a modifc


If we organized all variants by their locus, it is an easy 




Generating POS and CIGARS
-------------------------
    
    
    
    
## "Buffer bases" between simulated variants <a name="bufferbases"></a>

If we have variants adjacent to each other the most parsimonious description of the resulting variation can be different
from the original variants.

For example, considering `M=ATCGATCG` and an insertion and deletion as follows

    POS REF ALT
    1   A   ACC
    1   AT  A

We get

         1  2345678
    R    A  TCGATCG
    M    ACC CGATCG

This can, actually, be most parsimoniously expressed as

         1 2345678
    R    A TCGATCG
    M    ACCCGATCG

Which is a single base insertion followed by a SNP

    POS REF ALT
    1   A  AC
    2   T  C

For reasons of such ambiguity Mitty places a minimum 1 base "buffer" between variants, making the generated variant
identical to the most parsimonious description.
    