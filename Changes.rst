Planned changes
* change CLIs to allow parameter files to be fed in from std in (piped as a file)
* vcf-write should have option to remove unzipped .vcf file
* variants only reads - need to properly handle end of variants
* finish implementing ad hoc filtering - het/hom filtering in standard population model

Possible changes
* In read simulator/plugins change the 'S' for the sequences into 'obj' (like for variants)
* See if we can abstract any of the read plugin code into a base class
* Have a post filter in variant generation to do simple things like prune het/hom variants, prune variants too close
  together etc.

--- 1.18.0.dev0 ---

2015.10.06
* Enhancement: Full chain upto indel accuracy plot now works
* Enhancement: Ad hoc post filters implemented in standard population model.
het/hom filters still need to be implemented
* Bugfix: Now have a function return empty read array. This fixes an issue with read array concatenation: If we asked for
reads from variants only, but there were no variants, we would try to concatenate an empty list which would lead to
an error. This also fixes the problem that in such a condition the paired-endedness of the file would be uncertain.

2015.10.05
* Read length information added to qname (1.16.0.dev0)

2015.10.01
* Instead of making several different files write out the alignment accuracy in the original BAM itself.
  Still produce a perfect BAM as needed

2015.09.29
* Modified read simulator to allow reads to be generated over a sub-region of a chromosome.
  Coverage is correct. Sub-regions have to be set chromosome-by-chromosome.
  Parameter file format change is backwards compatible. Existing parameter files will work correctly with new version
* Added flag in read simulator to write gzipped fasta file.
  Existing parameter files will work correctly with new version