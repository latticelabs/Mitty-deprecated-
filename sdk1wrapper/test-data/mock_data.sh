#!/bin/bash
genomes generate -v -p variants.json
genomes genome-file write-vcf --sample-name 'g0_s0' genomes.h5 mock.vcf
reads generate -v -p reads.json
cp mock.* ../../../test-data/
cp chimera.fa.gz ../../../test-data/mock.fa.gz