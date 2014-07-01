"""Given the wrongly aligned reads and the reference genome plot the misalignment rate at each location. We really only
use this for haploid human data.

Usage:
  plot_bad_alignments.py --file=F

Options:
  --file=F       File with read analysis done
"""

import matplotlib
matplotlib.use('Agg')
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

if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__)

  data = numpy.load(args['--file'])
  ra = data['read data']
  hdr = data['reference data']
  chrom_keys_sorted = sorted(data['reference data'].keys(), key=lambda k: k[0])

  pylab.figure(figsize=(10, 1.25*len(chrom_keys_sorted)))
  pylab.subplots_adjust(left=0.1, right=0.99, top=0.99, bottom=0.01, hspace=0.02)
  for n, key in enumerate(chrom_keys_sorted):
    # For each chromosome, plot a histogram of incorrectly aligned reads by position
    pylab.subplot(len(chrom_keys_sorted), 1, n + 1)
    idx_wrong, = numpy.nonzero(~ra.correctly_aligned & (ra.correct_chrom_no == key[0]) & (ra.correct_chrom_copy == key[1]))
    num_total = numpy.count_nonzero((ra.correct_chrom_no == key[0]) & (ra.correct_chrom_copy == key[1]))
    seq_len = hdr[key]['len']
    if idx_wrong.size:
      pylab.hist(ra[idx_wrong].correct_pos, range=[0, seq_len], bins=100, normed=True, histtype='step')
    pylab.plot([centromere_pos[n+1], centromere_pos[n+1]], [0, pylab.getp(pylab.gca(), 'ylim')[1]], 'k:', lw=5)
    pylab.ylabel('{:s}\n({:2.2f}%)\n{:d}'.format(key, 100 * float(idx_wrong.size)/num_total, num_total))
    pylab.setp(pylab.gca(), xlim=[0, seq_len], xticks=[], yticks=[])

  fig_fname = os.path.splitext(args['--file'])[0] + '_analyzed.pdf'
  pylab.savefig(fig_fname)
  print 'Saved to {:s}'.format(fig_fname)
