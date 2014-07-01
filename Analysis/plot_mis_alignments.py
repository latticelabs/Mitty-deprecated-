"""Given the read analysis file, break the genome into bins and show a genome-wide histogram of where reads from each
bin have gone. Create as an animation, unless a single frame is specified.

Usage:
  plot_mis_alignments.py --file=F

Options:
  --file=F       File with read analysis done
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.animation as manimation
import numpy
import pylab
import docopt
import os


def compute_all(ref_data, read_data):
  idx, = numpy.nonzero(~read_data.correctly_aligned)
  correct_pos_whole_genome = read_data.correct_pos[idx].astype('u4')
  aligned_pos_whole_genome = read_data.aligned_pos[idx].astype('u4')

  offsets = [0] * 24
  for n in range(1, 24):
    offsets[n] = offsets[n - 1] + ref_data[(n, 1)]['len']

  for n in range(24):
    this_idx, = numpy.nonzero(read_data.correct_chrom_no[idx] == n + 1)
    correct_pos_whole_genome[this_idx] += offsets[n-1]

    this_idx, = numpy.nonzero(read_data.aligned_chrom_no[idx] == n + 1)
    aligned_pos_whole_genome[this_idx] += offsets[n-1]

  return correct_pos_whole_genome, aligned_pos_whole_genome


def plot_misaligned_hist(corr_pos, aligned_pos, genome_len, window_st=250e6, window_end=350e6, bins=50, xlim=None, ymax=3e-9):
  xlim = xlim or [0, genome_len]
  idx, = numpy.nonzero((corr_pos > window_st) & (corr_pos < window_end))
  pylab.hist(aligned_pos[idx], bins=bins, normed=True, range=xlim, histtype='step')
  pylab.plot([window_st, window_end], [0, 0], 'k-', lw=5)
  pylab.ylabel('Density')
  pylab.setp(pylab.gca(), xlim=xlim, ylim=[0, ymax], xticks=[], yticks=[])


# def plot_frame(corr_pos, aligned_pos, genome_len, window_st=250e6, window_end=350e6, bins=50, ymax=3e-9, ax=None):
#   pylab.axes(ax[0])
#   plot_misaligned_hist(corr_pos, aligned_pos, genome_len, window_st, window_end, bins, ymax=ymax)
#   pylab.axes(ax[1])
#   plot_misaligned_hist(corr_pos, aligned_pos, genome_len, window_st, window_end, bins, xlim=None, ymax=ymax)


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__)

  data = numpy.load(args['--file'])
  ra = data['read data']
  hdr = data['reference data']

  corr_pos, align_pos = compute_all(hdr, ra)

  FFMpegWriter = manimation.writers['ffmpeg']
  metadata = {'title': 'Misaligned reads', 'artist': 'Matplotlib', 'comment': 'Where do the misaligned reads go'}
  writer = FFMpegWriter(fps=15, metadata=metadata)

  fig, ax = pylab.subplots(1, 1, figsize=(14,4))
  pylab.subplots_adjust(left=0.05, right=0.99, bottom=0.13, top=0.95)
  pylab.xlabel('Genome pos')
  pylab.ylabel('Density')

  window_width = 100e6
  bins = 100
  with writer.saving(fig, "writer_test.mp4", 100):
    for read_st in range(0, int(3.14e9), int(window_width/10)):
      plot_misaligned_hist(corr_pos, align_pos, genome_len=3.14e9,
                           window_st=read_st, window_end=read_st+window_width, bins=bins, ymax=4e-9)
      writer.grab_frame()
      pylab.gca().cla()