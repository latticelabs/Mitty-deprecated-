"""In default mode, this script reads in a BAM file created by `reads.py` and aligns the reads by using the coordinates
stored in the seq id string. The whole genome file used to generate the reads OR the reference whole genome file is also
required. If the actual genome is given then exact alignments will be generated based on the actual POS data. If the
reference genome is given then all mappable reads will be aligned with the encoded POS and CIGAR strings. This aligned
bam file can be read in using a visualizer to debug the simulation chain.

If we use the 'analyze' option then the input BAM file is treated as the output of an aligner and the reads are organized
into a table that allows us to analyze the performance of the aligner.

Usage:
cheata  --inbam=INBAM  --outbam=OUTBAM  --wg=WG  [-v]
cheata  analyze --inbam=INBAM  --wg=WG [-v]

Options:
  analyze                 Instead of aligning the BAM file, analyze its alignment
  --inbam=INBAM           Input bam file name of reads
  --outbam=OUTBAM         Output (perfectly aligned) bam file name
  --wg=WG                 Whole genome file. Needed for the index, to match BAM seq ids to chrom_key in read qname
  -v                      Dump detailed logger messages

Notes:
1. Recall that the seq id of each read is the string 'ch:cp|rN|POS1|CIGAR1|PABS1|POS2|CIGAR2|PABS2'
   Unpaired reads have no POS2, CIGAR2
2. As of v 0.2.0 there is a quirk with pysam writing BAM files where alignments have some hidden problem when written as
   BAM from cheata. My work around is to save as SAM and then convert to BAM which seems to not create any problems.
3. The `split` mode enables us to break the larger aligned .bam file into smaller files. Since we are often interested
   in just the errors, we can look at just the _wrong.bam which is often much smaller. We also get a summary of the
   alignment performance for easy plotting.
"""
__version__ = '0.2.0'

import os  # Needed for filesize (for sequence length) and messing with filenames
import pysam  # Needed to read/write BAM files
import numpy
import cPickle
import docopt
import h5py
from reads import interpret_read_qname
import logging
logger = logging.getLogger(__name__)


def sort_and_index_bam(bamfile):
  """Do the filename gymnastics required to end up with a sorted, indexed, bam file."""
  # samtools sort adds a '.bam' to the end of the file name.
  pysam.sort(bamfile, os.path.splitext(bamfile)[0])
  pysam.index(bamfile)


def align(inbam, outbam, ref_seqs=[], absolute=True):
  """ref_seqs - [(seq_id, seq_len) ...]  in correct order so they can be recovered from correct_chrom_no"""
  in_bamfile = pysam.Samfile(inbam, 'rb')
  out_hdr = in_bamfile.header
  out_hdr['SQ'] = [{'SN': seq_id.split(' ')[0].strip(), 'LN': seq_len} for (seq_id, seq_len) in ref_seqs]
  out_bamfile = pysam.Samfile(outbam, 'wb', header=out_hdr)

  for cnt, read in enumerate(in_bamfile):
    correct_chrom_no, correct_chrom_copy, correct_pos, correct_cigar, absolute_pos = interpret_read_qname(read)
    read.mapq = 100  # It's better to set this
    read.tid = correct_chrom_no - 1
    if absolute:
      read.pos = absolute_pos
      read.cigarstring = '{:d}M'.format(len(read.sequence))
    else:
      read.pos = correct_pos  # 0-indexed
      read.cigarstring = correct_cigar
    out_bamfile.write(read)

    if cnt % 100000:
      logger.debug('{:d} reads done'.format(cnt))

  logger.debug('{:d} reads done'.format(cnt))
  out_bamfile.close()
  sort_and_index_bam(outbam)


def analyze(inbam, ref_seqs=[], data_save_fp):
  """Split the aligned BAM file into three sets of reads - properly aligned, misaligned and unmapped - and save them in
  separate bam files
  """
  logger.error('Not updated yet')
  raise RuntimeError

  #The rname in the aligned read corresponds to the sequence in the list in_bam.header['SQ']
  #We need to match that up with the chromosome code we have in the wg file
  seq_name_list = in_bam.references

  total_reads = in_bam.mapped + in_bam.unmapped   #total_reads = int(pysam.view('-c', in_bam.filename)[0].strip())
  logger.debug('Total reads: {:d}'.format(total_reads))

  #This array stores our parsed read data
  # read_analysis = numpy.recarray(shape=(total_reads,),
  #                                dtype=[('correct_chrom', '5S'), ('correct_pos', '<i4'),
  #                                       ('correct_cigar', '255S'),
  #                                       ('aligned_chrom', '5S'), ('aligned_pos', '<i4'),
  #                                       ('aligned_cigar', '255S'),
  #                                       ('mapped', bool), ('correctly_aligned', bool),
  #                                       ('mapping_qal', '<i1')])

  read_analysis = numpy.recarray(shape=(total_reads,),
                                 dtype=[('correct_chrom_no', '<i1'), ('correct_chrom_copy', '<i1'), ('correct_pos', '<u4'),
                                        ('aligned_chrom_no', '<i1'), ('aligned_chrom_copy', '<i1'), ('aligned_pos', '<u4'),
                                        ('absolute_pos', 'u4'),
                                        ('mapped', bool), ('correctly_aligned', bool),
                                        ('mapping_qual', '<i1')])

  for n, read in enumerate(in_bam):
    correct_chrom_no, correct_chrom_copy, correct_pos, correct_cigar, absolute_pos = interpret_read_qname(read)#reads.interpret_read_qname(read)
    read_analysis[n].correct_chrom_no = correct_chrom_no
    read_analysis[n].correct_chrom_copy = correct_chrom_copy
    read_analysis[n].correct_pos = correct_pos
    #read_analysis[n].correct_cigar = correct_cigar
    if read.flag & 0x4:  # Unmapped
      read_analysis[n].mapped = False
      unmapped_bam.write(read)
    else:
      read_analysis[n].mapped = True
      aligned_chrom_no, aligned_chrom_copy = rev_idx[seq_name_list[read.rname]]
      aligned_pos = read.pos + 1  # PySam uses 0 indexing ...

      read_analysis[n].aligned_chrom_no = aligned_chrom_no
      read_analysis[n].aligned_chrom_copy = aligned_chrom_copy
      read_analysis[n].aligned_pos = aligned_pos
      #read_analysis[n].correct_cigar = correct_cigar
      read_analysis[n].mapping_qual = read.mapq

      if correct_chrom_no == aligned_chrom_no and correct_chrom_copy == aligned_chrom_copy and correct_pos == aligned_pos:
        correct_bam.write(read)
        read_analysis[n].correctly_aligned = True
      else:
        wrong_bam.write(read)
        read_analysis[n].correctly_aligned = False

  cPickle.dump({
    'read data': read_analysis,
    'reference data': {rev_idx[seq_n]: {'name': seq_n, 'len': l} for seq_n, l in zip(seq_name_list, in_bam.lengths)}
  }, data_save_fp, protocol=cPickle.HIGHEST_PROTOCOL)

  logger.debug('{:d} reads done'.format(total_reads))
  good_cnt = numpy.count_nonzero(read_analysis.correctly_aligned)
  bad_cnt = in_bam.mapped - good_cnt
  logger.debug('Found {:d} good reads, {:d} bad reads, {:d} unmapped reads'.format(good_cnt, bad_cnt, in_bam.unmapped))


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  with h5py.File(args['--wg'], 'r') as ref_fp:
    align(inbam, outbam, ref_seqs=[], absolute=True)

  if args['split']:  # We are in split_good_bad mode
    inbam_name = args['--inbam']
    correct_bam_name = os.path.splitext(inbam_name)[0] + '_correct.bam'
    wrong_bam_name = os.path.splitext(inbam_name)[0] + '_wrong.bam'
    unmapped_bam_name = os.path.splitext(inbam_name)[0] + '_unmapped.bam'
    data_file_name = os.path.splitext(inbam_name)[0] + '_alignment_data.pkl'
    with genome.WholeGenome(fname=args['--wg']) as wg:
      in_bam_f=pysam.Samfile(inbam_name, 'rb')
      split(wg, in_bam=in_bam_f,
            correct_bam=pysam.Samfile(correct_bam_name, 'wb', header=in_bam_f.header),
            wrong_bam=pysam.Samfile(wrong_bam_name, 'wb', header=in_bam_f.header),
            unmapped_bam=pysam.Samfile(unmapped_bam_name, 'wb', header=in_bam_f.header),
            data_save_fp=open(data_file_name, 'wb'))
      logger.debug('Saved data to:\n{:s}\n{:s}\n{:s}\n{:s}'.format(correct_bam_name, wrong_bam_name, unmapped_bam_name, data_file_name))
  else:  # We are in align mode



    seq_name, seq_len = open(args['--heada']).readlines()
    align(args['--inbam'], args['--outbam'], seq_name, seq_len)

