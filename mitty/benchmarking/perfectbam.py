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

import os
import time
import io

import pysam
from collections import Counter
import click

import mitty.benchmarking.creed as creed
import mitty.lib.io as mio  # For the bam sort and index function

import logging
logger = logging.getLogger(__name__)


def main(bam_in_fp, bam_out_fp=None, db_name=None, cr_fname=None, window=0, extended=False, show_progress_bar=False):
  """Main processing function that goes through the bam file, analyzing read alignment and writing out

  :param bam_in_fp:  Pointer to original BAM
  :param bam_out_fp: Pointer to perfect BAM being created
  :param db_name:    database name
  :param cr_fname:   name of file to store analysed reads in
  :param window:     Tolerance window for deciding if read is correctly aligned
  :param extended:   If True write out new style CIGARs (With '=' and 'X')
  :param show_progress_bar: Our famous progress bar, if we want it
  :return: number of reads processed
  """
  debug_db = creed.ReadDebugDB(db_name) if db_name else None
  cat_reads = creed.CategorizedReads(fname=cr_fname, seq_names=bam_in_fp.references, seq_lengths=bam_in_fp.lengths, copies=2) if cr_fname else None

  total_read_count = bam_in_fp.mapped + bam_in_fp.unmapped
  n0 = 0
  delta_read_cnt = int(0.01 * total_read_count)

  total_reads_cntr, incorrectly_aligned_reads_cntr, unmapped_reads_cntr = Counter(), Counter(), Counter()

  analyze_read = creed.analyze_read
  cra = cat_reads.append if cat_reads else lambda **kwargs: None  # The lambda is a NOP

  with click.progressbar(length=total_read_count, label='Processing BAM', file=None if show_progress_bar else io.BytesIO()) as bar:
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

      cra(chrom=chrom, cpy=cpy, pos=pos, read_len=read.query_length, code=error_type)

      if bam_out_fp:
        # Now write out the perfect alignment
        read.is_reverse = 1 - ro
        read.mate_is_reverse = ro
        read.mate_is_unmapped = False  # Gotta check this - what if mate is deep in an insert?
        read.reference_id = chrom - 1
        read.pos = pos
        read.cigarstring = cigar  # What if this is deep in an insert?
        bam_out_fp.write(read)

      if n - n0 > delta_read_cnt:
        bar.update(n - n0)
        n0 = n

  if debug_db:
    debug_db.write_summary_to_db(total_reads_cntr, incorrectly_aligned_reads_cntr, unmapped_reads_cntr, bam_in_fp.header['SQ'])
    debug_db.commit_and_create_db_indexes()

  if cat_reads: cat_reads.finalize()

  return int(total_read_count)

@click.command()
@click.version_option()
@click.argument('inbam', type=click.Path(exists=True))
@click.option('--perfectbam', help='Write out perfect BAM to this file', type=click.Path())
@click.option('--debugdb', help='Write mis aligned reads to this database', type=click.Path())
@click.option('--catreads', help='Write categorized reads to this file', type=click.Path())
@click.option('--window', help='Size of tolerance window', default=0, type=int)
@click.option('-x', is_flag=True, help='Use extended CIGAR ("X"s and "="s) rather than traditional CIGAR (just "M"s)')
@click.option('-v', count=True, help='Verbosity level')
@click.option('-p', is_flag=True, help='Show progress bar')
# [--perfectbam=PBAM] [--debugdb=DB] [--catreads=CR] [--window=WN] [-x] [-v] [-p]
def cli(inbam, perfectbam, debugdb, catreads, window, x, v, p):
  """Command line script entry point."""
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  if perfectbam is None and debugdb is None and catreads is None:
    print('No outputs specified. Easiest gig ever.')
    return

  bam_in_fp = pysam.AlignmentFile(inbam, 'rb')
  bam_out_fp = pysam.AlignmentFile(perfectbam, 'wb', template=bam_in_fp) if perfectbam else None

  if debugdb:
    try:
      os.remove(debugdb)
    except OSError:
      pass

  if catreads:
    try:
      os.remove(catreads)
    except OSError:
      pass

  t0 = time.time()
  read_count = main(bam_in_fp=bam_in_fp, bam_out_fp=bam_out_fp, db_name=debugdb, cr_fname=catreads,
                    window=window, extended=x, show_progress_bar=p)
  if bam_out_fp: bam_out_fp.close()
  t1 = time.time()
  logger.debug('Analyzed {:d} reads in BAM in {:2.2f}s'.format(read_count, t1 - t0))

  if perfectbam:
    t0 = time.time()
    mio.sort_and_index_bam(perfectbam)
    t1 = time.time()
    logger.debug('Sort and indexed perfect BAM in {:2.2f}s'.format(t1 - t0))


if __name__ == "__main__":
  cli()