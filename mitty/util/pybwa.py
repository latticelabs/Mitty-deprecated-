"""A simple wrapper round BWA because I got tired of writing out the invocations

  bwa mem -t 8 -p ../data/red_alga.fa.gz  Out/reads_c.fq > Out/temp.sam
  samtools view -Sb  Out/temp.sam > Out/temp.bam
  samtools sort Out/temp.bam  Out/reads
  samtools index Out/reads.bam

Usage:
  pybwa <ref> <fastq> <bam> [-S] [-P] [-t=T] [-p] [-v...]

Options:
  <ref>     fasta.gz reference
  <fastq>   fastq file
  <bam>     BAM file to write to
  -t=T      Number of threads [default: 8]
  -p        interleaved paired end file?
  -S        skip mate rescue
  -P        skip pairing
  -v        Verbose: print commands and messages from commands. Otherwise, suppress ALL output incl. stderr
            The more -v s the more verbose
"""
import subprocess
import tempfile
import os

import docopt


def cli():
  """Script entry point"""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__)
  _, sam_file = tempfile.mkstemp('.sam')
  _, temp_bam_file = tempfile.mkstemp('.bam')
  verb = args['-v']

  commands = [
    'bwa mem -t {t} {p} {S} {P} {ref} {fastq} > {sam}'.
      format(t=int(args['-t']), p='-p' if args['-p'] else '', S='-S' if args['-S'] else '', P='-P' if args['-P'] else '',
             ref=args['<ref>'], fastq=args['<fastq>'], sam=sam_file),
    'samtools view -Sb {:s} > {:s}'.format(sam_file, temp_bam_file),
    'samtools sort {:s} -o {:s}'.format(temp_bam_file, args['<bam>']),
    'samtools index {:s}'.format(args['<bam>']),
  ]
  for cmd in commands:
    if verb: print(cmd)
    if subprocess.call(cmd + (' 2>/dev/null' if verb < 2 else ''), shell=True):
      print('Some error, quitting')
      break

  os.remove(sam_file)
  os.remove(temp_bam_file)

if __name__ == '__main__':
  cli()