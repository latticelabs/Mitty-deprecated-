Planned changes
* vcf-write should have option to remove unzipped .vcf file
* variants only reads - need to properly handle end of variants

2015.10.01
* Instead of making several different files write out the alignment accuracy in the original BAM itself.
  Still produce a perfect BAM as needed

2015.09.29
* Modified read simulator to allow reads to be generated over a sub-region of a chromosome.
  Coverage is correct. Sub-regions have to be set chromosome-by-chromosome.
  Parameter file format change is backwards compatible. Existing parameter files will work correctly with new version
* Added flag in read simulator to write gzipped fasta file.
  Existing parameter files will work correctly with new version