#!/bin/bash
set -x
: Create whole genome file from fasta files
python fasta2wg.py --index=Examples/Genome/wg_chimera.json --wg=Examples/Genome/Out/chimera.h5 --fa=Examples/Genome/Out/chimera.fa.gz -v
: Use the describe command to reveal the contents of the genome file
python fasta2wg.py describe --wg=Examples/Genome/Out/chimera.h5