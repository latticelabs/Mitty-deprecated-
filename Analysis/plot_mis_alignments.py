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

centromere_pos = {
  1: 125.0e6,
  2: 93.3e6,
  3: 91.0e6,
  4: 50.4e6,
  5: 48.4e6,
  6: 61.0e6,
  7: 59.9e6,
  8: 45.6e6,
  9: 49.0e6,
  10: 40.2e6,
  11: 53.7e6,
  12: 35.8e6,
  13: 17.9e6,
  14: 17.6e6,
  15: 19.0e6,
  16: 36.6e6,
  17: 24.0e6,
  18: 17.2e6,
  19:	26.5e6,
  20: 27.5e6,
  21: 13.2e6,
  22: 14.7e6,
  23: 60.6e6,
  24: 12.5e6
}


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

  return correct_pos_whole_genome, aligned_pos_whole_genome, offsets


def plot_misaligned_hist(corr_pos, aligned_pos, chrom_pos, genome_len, window_st=250e6, window_end=350e6, bins=50, xlim=None, ymax=3e-9):
  xlim = xlim or [0, genome_len]
  idx, = numpy.nonzero((corr_pos > window_st) & (corr_pos < window_end))
  pylab.hist(aligned_pos[idx], bins=bins, normed=False, range=xlim, histtype='step')
  for p in chrom_pos:
    pylab.plot([p, p], [0, ymax], 'k:', lw=2)
  pylab.plot([window_st, window_end], [0, 0], 'k-', lw=5)
  pylab.xlabel('Genome pos')
  pylab.ylabel('Density\n({:d} misalinged reads)'.format(idx.size))
  pylab.setp(pylab.gca(), xlim=xlim, ylim=[0, ymax], xticks=[], yticks=[])


def plot_misalignment_hexbin(corr_pos, aligned_pos, window_st=250e6, window_end=350e6, bins=50):
  window_w = window_end - window_st
  x0 = window_st - 5 * window_w
  x1 = window_end + 5 * window_w
  idx, = numpy.nonzero((corr_pos > x0) & (corr_pos < x1) & (aligned_pos > x0) & (aligned_pos < x1))
  pylab.hexbin(corr_pos[idx], aligned_pos[idx], gridsize=bins, cmap=pylab.cm.gray_r, extent=(x0, x1, x0, x1))
  pylab.plot((x0, x1), (x0, x1), 'k:')
  pylab.axis('scaled')
  pylab.setp(pylab.gca(), xlim=(x0, x1), ylim=(x0, x1))


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

  corr_pos, align_pos, chrom_pos = compute_all(hdr, ra)

  FFMpegWriter = manimation.writers['ffmpeg']
  metadata = {'title': 'Misaligned reads', 'artist': 'Matplotlib', 'comment': 'Where do the misaligned reads go'}
  writer = FFMpegWriter(fps=10, metadata=metadata)

  fig, ax = pylab.subplots(2, 1, figsize=(14, 9))
  pylab.subplots_adjust(left=0.05, right=0.99, bottom=0.13, top=0.95)

  window_width = 50e6
  window_advance = int(window_width / 20.)
  bins = 100
  hexbins = 50
  ymax = 10000#1e-8
  genome_len = 3.08e9
  movie_to = int(genome_len)
  with writer.saving(fig, "writer_test.mp4", 150):
    for read_st in range(0, int(movie_to), window_advance):
      pylab.axes(ax[0])
      plot_misaligned_hist(corr_pos, align_pos, chrom_pos=chrom_pos, genome_len=genome_len,
                           window_st=read_st, window_end=read_st+window_width, bins=bins, ymax=ymax)
      pylab.axes(ax[1])
      plot_misalignment_hexbin(corr_pos, align_pos, window_st=read_st, window_end=read_st+window_width, bins=hexbins)
      writer.grab_frame()
      ax[0].cla()
      ax[1].cla()