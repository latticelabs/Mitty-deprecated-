"""Plot mis-alignments as circular plots or 2D histograms

Commandline::

  Usage:
    plot_gc_bias <fa.gz> <bam> [--out=OUT] [--min_map_qual=QUAL] [--win=WIN] [--g0=G0] [--g1=G1] [pdf] [-p]

  Options:
    <fa.gz>    Name of fasta (reference) file
    <bam>      Name of BAM file
    --out=OUT  Output file name
    --min_map_qual=QUAL  Minimum mapping quality of read to accept [default: 20]
    --win=WIN  gc computation window [default: 10000]
    --g0=G0    lower gc limit of plot [default: 0.3]
    --g1=G1    upper gc limit of plot [default: 0.7]
    pdf        Should the output be pdf? (png otherwise) (Same as ending file name in --out with .pdf)
    -p         Display progress bar
"""

import docopt
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pysam

from mitty.lib import progress_bar
import mitty.lib.io as mio


def gc_content_of_genome(ref, gc_window=10000):
  """Compute GC content of all sequences in bins

  :param ref:       Fasta object (should be loaded on init)
  :param gc_window: length of window for gc computation
  :return: list of np arrays
  """
  gc = []
  for r in ref:
    gc += [gc_content(r['seq'], gc_window)]
  return gc


def gc_content(seq, gc_window=10000):
  """Compute GC content of sequence

  :param seq: sequence
  :param gc_window: length of window for gc computation
  :return: np array
  """
  f_gc_window = float(gc_window)
  n_bins = int(np.ceil(len(seq) / f_gc_window))
  this_gc = np.empty(n_bins, dtype=float)
  for n in range(n_bins):
    st = n * gc_window
    stp = st + gc_window
    this_gc[n] = (seq.count('G', st, stp) + seq.count('C', st, stp)) / f_gc_window
  return this_gc


def process_bam_for_gc(bam_in_fp, ref, gc_window=10000, min_map_qual=20, progress_bar_func=None):
  """Goes through the bam file, analyzing gc content of reads.

  :param bam_in_fp: pointer to opened BAM file
  :param ref: Fasta object
  :param gc_window: The size of the gc_window used
  :param progress_bar_func: if a proper progress bar function is passed, this will show a progress bar as we complete
  :return: x, y two arrays with gc content ratio (x) and read counts (y) for those bins
  """
  #gc_counts = [np.zeros(s.shape[0], dtype=int) for s in gc]
  gc = {}
  gc_counts = {}
  total_mapped_read_count = float(bam_in_fp.mapped)
  reads_processed = 0
  f0 = 0
  if progress_bar_func is not None: progress_bar_func('Processing BAM ', 0, 80)
  for read in bam_in_fp:
    if read.is_unmapped: continue
    if read.mapping_quality < min_map_qual: continue
    if read.reference_id not in gc_counts:
      gc[read.reference_id] = gc_content(ref[read.reference_id + 1]['seq'], gc_window)
      gc_counts[read.reference_id] = np.zeros(gc[read.reference_id].shape[0], dtype=int)
    gc_counts[read.reference_id][int(read.pos / gc_window)] += 1
    reads_processed += 1
    if progress_bar_func is not None:
      f = reads_processed / total_mapped_read_count
      if f - f0 >= 0.01:
        progress_bar_func('Processing BAM ', f, 80)
        f0 = f
  if progress_bar_func is not None: print('\n')

  return gc, gc_counts


def plot_gc_bias(gc_list, gc_counts_list, xlim):
  for k in gc_counts_list.keys():
    plt.scatter(gc_list[k], gc_counts_list[k], s=1, marker='.')
  plt.setp(plt.gca(), xlim=xlim, xlabel='GC ratio', ylabel='Read count')


def cli():
  """Command line script entry point."""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__)

  ref = mio.Fasta(multi_fasta=args['<fa.gz>'])
  with pysam.AlignmentFile(args['<bam>'], 'rb') as bam_in_fp:
    #t0 = time.time()
    #gc = gc_content_of_genome(ref, gc_window=int(args['--win']))
    gc, gc_counts = process_bam_for_gc(bam_in_fp, ref, gc_window=int(args['--win']),
                                       min_map_qual=int(args['--min_map_qual']),
                                       progress_bar_func=progress_bar if args['-p'] else None)
    plot_gc_bias(gc, gc_counts, xlim=[float(args['--g0']), float(args['--g1'])])
    plot_suffix = '.pdf' if args['pdf'] else '.png'
    fname = args['--out'] or args['<bam>'] + '_gc_bias' + plot_suffix
    plt.savefig(fname)

if __name__ == '__main__':
  cli()