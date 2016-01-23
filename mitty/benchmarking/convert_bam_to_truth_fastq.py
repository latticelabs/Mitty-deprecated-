import io
import time

import click
import pysam

import logging
logger = logging.getLogger(__name__)


@click.command()
@click.argument('inbam', type=click.Path(exists=True))
@click.argument('fastq', type=click.Path())
@click.option('--mq-threshold', type=click.IntRange(min=0, max=255), help='Skip reads below this threshold', default=0)
@click.option('-v', count=True, help='Verbosity level')
@click.option('-p', is_flag=True, help='Show progress bar')
def cli(inbam, fastq, mq_threshold, v, p):
  """This consumes a BAM, treats it as a truth dataset and writes out all the mapped reads."""
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  bam_in_fp = pysam.AlignmentFile(inbam, 'rb')
  total_read_count = bam_in_fp.mapped + bam_in_fp.unmapped  # Sadly, this is only approximate
  progress_bar_update_interval = int(0.01 * total_read_count)

  cnt, rs = 0, 0
  t0 = time.time()
  with click.progressbar(length=total_read_count, label='Processing BAM',
                         file=None if p else io.BytesIO()) as bar, open(fastq, 'w') as fastq_out_fp:
    for cnt, rs in process_file(bam_in_fp=bam_in_fp, fastq_out_fp=fastq_out_fp,
                                mq_threshold=mq_threshold,
                                progress_bar_update_interval=progress_bar_update_interval):
      bar.update(progress_bar_update_interval)
  t1 = time.time()
  logger.debug('Analyzed {:d} reads in {:2.2f}s. Wrote out {:d} templates'.format(cnt, t1 - t0, rs))



def process_file(bam_in_fp, fastq_out_fp, mq_threshold, progress_bar_update_interval=100):
  """Main processing function that goes through the bam file, analyzing read alignment and writing out

  :param bam_in_fp:  Pointer to original BAM
  :param fastq_out_fp: Pointer to file to be written
  :param mq_threshold: Pointer to PERBAM being created
  :param progress_bar_update_interval: how many reads to process before yielding (to update progress bar as needed)
  :return: number of reads processed
  """
  n0 = progress_bar_update_interval
  read_cache = {}  # dictionary with qname as key and read as value
  read_serial = 0
  for tot_read_cnt, read in enumerate(bam_in_fp):
    n0 -= 1
    if n0 == 0:
      yield tot_read_cnt, read_serial
      n0 = progress_bar_update_interval

    if read.is_unmapped or read.mate_is_unmapped:
      continue

    if read.reference_id != read.next_reference_id:  # Mates are in different chroms, don't want that
      continue

    if read.is_paired:
      if read.qname in read_cache:  # Yay we found the mate
        read2 = read_cache.pop(read.qname)
        if read.mapping_quality >= mq_threshold and read2.mapping_quality >= mq_threshold:
          read_serial = flush_reads([read, read2], fastq_out_fp, read_serial)
      else:
        read_cache[read.qname] = read
    else:
      if read.mapping_quality >= mq_threshold:
        read_serial = flush_reads([read], fastq_out_fp, read_serial)

  yield tot_read_cnt + 1, read_serial  # tot_read_cnt starts from 0 actually ...


def flush_reads(read_l, fastq_out_fp, read_serial):
  """Write out reads to a fastq file.
  qname is 'read_serial|chrom|copy|ro|pos|rlen|cigar|ro|pos|rlen|cigar'

  :param read_l:
  :param fastq_out_fp:
  :param read_serial:
  :return: read_serial, updated
  """
  qname = '|'.join([str(read_serial), str(read_l[0].reference_id + 1), '0'])
  for ro, r in enumerate(read_l):
    qname += '|'.join([str(ro), str(r.pos), str(r.query_length), r.cigarstring])

  for n, r in enumerate(read_l):
    fastq_out_fp.write('@' + qname + ('/' + str(n + 1) if len(read_l) > 0 else '') + '\n'
                       + r.query_sequence + '\n+\n'
                       + '~' * len(r.query_sequence) + '\n')

  return read_serial + 1