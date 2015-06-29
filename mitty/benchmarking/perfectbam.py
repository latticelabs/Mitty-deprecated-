"""In order to perform our alignment analysis we need to go through the BAM file and classify whether each read was
correctly aligned or not. The following back-of-the-envelope calculations indicate to us what size of data we are
expecting for a 30x illumina whole genome dataset:

# reads = 30 * 3e9 / 100 = 900e6 (900 million reads)

If about 5% if these are misaligned/unmapped, we have 45 million misaligned reads

For debugging purposes we store complete read data for only the incorrectly aligned reads. These are stored in a
database which offers efficient enough retrieval for debugging.

For comprehensive data analysis, where we want to compute alignment accuracy over different spots of the genome, but are
not interested in the read details - like sequence and base quality score - we only store the following information

chrom,
cpy,
start_pos  - 4 bytes  start of read
stop_pos   - 4 bytes  end of read    = start of read + read length
category   - 1 byte
---------------------
             9 bytes

The category flags are defined as follows

---------------------------------
| 7 | 6 | 5 | 4 | 3 | 2 | 1 | 0 |
---------------------------------
          |   |   |   |   |   |
          |   |   |   |   |   \------  1 => chrom was wrong
          |   |   |   |   \----------  1 => pos was wrong
          |   |   |   \--------------  1 => cigar was wrong
          |   |   \------------------  1 => unmapped
          |   \----------------------  1 => read from reference region (no variants)
          \--------------------------  1 => mate is from reference region


To store this data for a 30x WG Illumina run will take 9 x 900e6 ~ 8.1 GB

This can  be loaded onto a modern machine. For this reason, we don't do anything very fancy and simply store the
categorized reads as a set of compressed numpy arrays using readily available mechanisms.

The categorized reads are stored in POS sorted order. The primary analysis we are are considering is to count the number
of correct and total reads in the vicinity of a series of positions (usually variants) or for the whole sequence without
regard to position (e.g. for reference reads).

The most efficient way to do this is to walk through the sorted list of reads side by side with the sorted list of
variants and count the relevant reads as we go.
"""
from mitty.version import __version__
__cmd__ = """perfectbam ({:s}): Categorize read alignments in a BAM file

Commandline::

  Usage:
    perfectbam <inbam> [--perfectbam=PBAM] [--debugdb=DB] [--catreads=CR] [--window=WN] [-x] [-v] [-p]

  Options:
    <inbam>             Input bam file name of reads
    --perfectbam=PBAM   Perfect BAM will be written to this file
    --debugdb=DB        Detailed data about misaligned reads will be written to this file
    --catreads=CR       All categorized reads will be saved to this file
    --window=WN         Size of tolerance window [default: 0]
    -x                  Use extended CIGAR ('X's and '='s) rather than traditional CIGAR (just 'M's)
    -v                  Dump detailed logger messages
    -p                  Show progress bar
""".format(__version__)

__param__ = """Given a bam file containing simulated reads aligned by a tool
  1. Produce a new bam that re-aligns all reads so that their alignment is perfect

  2. Produce a database file containing data about the misaligned and unmapped reads with the following tables::

  summary::

      chrom           - chromosome number
      total_reads     - total reads on this chrom
      incorrect_reads - total incorrect reads on this chrom
      unmapped_reads  - unmapped reads
      seq_len         - len of this sequence
      seq_id          - full seq name as found in BAM header

  reads::

       qname              -  qname of the read
       error_type         -  type of error 3 bit number  bit 0=chrom, 1=pos, 2=cigar, 3=unmapped
       correct_chrom      -  correct chromosome number of read
       correct_pos        -  correct position of read
       correct_cigar      -  correct cigar of read
       aligned_chrom      -  actual aligned chromosome
       aligned_pos        -  actual aligned pos
       aligned_cigar      -  actual aligned cigar
       mapping_qual       -  mapping quality
       mate_is_unmapped   -  is mate unmapped
       seq                -  actual sequence string

  3. Produce a compressed numpy array file with details of all the reads in the file

      chrom_cc   - 8 bit int: bit 0 -> chrom copy 0=0, 1=1  bit 1 onwards -> original chromosome of read
      pos        - 32 bit int: original position of read
      code       - reference_read_bit | error_type
                   bit 4 = 1 if read is from reference

      We could have split up the arrays into no_of_chroms x 2 arrays and dropped the chromosome field for space savings
      but this makes the code more complex since we no longer no how large to make the arrays.
      It also makes computations on the whole genome faster.
  """
__doc__ = __cmd__ + '\n\nDetails:\n\n' + __doc__

import os
import time

import pysam
from collections import Counter
import docopt

import mitty.benchmarking.creed as creed
import mitty.lib.io as mio  # For the bam sort and index function
from mitty.lib import progress_bar

import logging
logger = logging.getLogger(__name__)


def main(bam_in_fp, bam_out_fp=None, db_name=None, cr_fname=None, window=0, extended=False, progress_bar_func=None):
  """Main processing function that goes through the bam file, analyzing read alignment and writing out

  :param bam_in_fp:  Pointer to original BAM
  :param bam_out_fp: Pointer to perfect BAM being created
  :param db_name:    database name
  :param cr_fname:   name of file to store analysed reads in
  :param window:     Tolerance window for deciding if read is correctly aligned
  :param extended:   If True write out new style CIGARs (With '=' and 'X')
  :param progress_bar_func: Our famous progress bar, if we want it
  :return: number of reads processed
  """
  debug_db = creed.ReadDebugDB(db_name) if db_name else None
  cat_reads = creed.CategorizedReads(fname=cr_fname, seq_names=bam_in_fp.references, seq_lengths=bam_in_fp.lengths, copies=2) if cr_fname else None

  total_read_count = float(bam_in_fp.mapped + bam_in_fp.unmapped)
  f0 = 0
  total_reads_cntr, incorrectly_aligned_reads_cntr, unmapped_reads_cntr = Counter(), Counter(), Counter()

  if progress_bar_func is not None: progress_bar_func('Processing BAM ', 0, 80)

  analyze_read = creed.analyze_read
  cra = cat_reads.append if cat_reads else lambda **kwargs: None  # The lambda is a NOP

  for n, read in enumerate(bam_in_fp):
    read_serial, chrom, cpy, ro, pos, cigar, error_type = analyze_read(read, window, extended)
    if read_serial is None: continue
    total_reads_cntr[chrom] += 1

    if error_type & 0b1111:
      if error_type & 0b1000:  # Unmapped read
        unmapped_reads_cntr[chrom] += 1
      else:
        incorrectly_aligned_reads_cntr[chrom] += 1
      if debug_db:
        debug_db.write_read_to_misaligned_read_table(
          [read_serial, read.qname, error_type, chrom, pos, cigar,
           read.reference_id + 1, read.pos, read.cigarstring,
           read.mapq, read.mate_is_unmapped, read.query_sequence])

    cra(chrom, cpy, pos, read.query_length, error_type)

    if bam_out_fp:
      # Now write out the perfect alignment
      read.is_reverse = 1 - ro
      read.mate_is_reverse = ro
      read.mate_is_unmapped = False  # Gotta check this - what if mate is deep in an insert?
      read.reference_id = chrom - 1
      read.pos = pos
      read.cigarstring = cigar  # What if this is deep in an insert?
      bam_out_fp.write(read)

    if progress_bar_func is not None:
      f = n / total_read_count
      if f - f0 >= 0.01:
        progress_bar_func('Processing BAM ', f, 80)
        f0 = f
  if progress_bar_func is not None: print('\n')

  if debug_db:
    debug_db.write_summary_to_db(total_reads_cntr, incorrectly_aligned_reads_cntr, unmapped_reads_cntr, bam_in_fp.header['SQ'])
    debug_db.commit_and_create_db_indexes()

  if cat_reads: cat_reads.finalize()

  return int(total_read_count)


def cli():
  """Command line script entry point."""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__cmd__, ['-h'])
  else:
    args = docopt.docopt(__cmd__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  if args['--perfectbam'] is None and args['--debugdb'] is None and args['--catreads'] is None:
    print('No outputs specified. Easiest gig ever.')
    return

  bam_in_fp = pysam.AlignmentFile(args['<inbam>'], 'rb')
  pbam_fname = args['--perfectbam']
  bam_out_fp = pysam.AlignmentFile(pbam_fname, 'wb', template=bam_in_fp) if pbam_fname else None
  db_name = args['--debugdb']
  cr_fname = args['--catreads']

  if db_name:
    try:
      os.remove(db_name)
    except OSError:
      pass

  if cr_fname:
    try:
      os.remove(cr_fname)
    except OSError:
      pass

  t0 = time.time()
  read_count = main(bam_in_fp=bam_in_fp, bam_out_fp=bam_out_fp, db_name=db_name, cr_fname=cr_fname,
                    window=int(args['--window']), extended=bool(args['-x']), progress_bar_func=progress_bar if args['-p'] else None)
  if bam_out_fp: bam_out_fp.close()
  t1 = time.time()
  logger.debug('Analyzed {:d} reads in BAM in {:2.2f}s'.format(read_count, t1 - t0))

  if pbam_fname:
    t0 = time.time()
    mio.sort_and_index_bam(args['--perfectbam'])
    t1 = time.time()
    logger.debug('Sort and indexed perfect BAM in {:2.2f}s'.format(t1 - t0))


if __name__ == "__main__":
  cli()