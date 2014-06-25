"""Given the wrongly aligned reads and the reference genome plot the misalignment rate at each location.

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

if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__)

  data = numpy.load(args['--file'])
  ra = data['read data']
  hdr = data['reference data']
  chrom_keys_sorted = sorted(data['reference data'].keys(), key=lambda k: int(k.split(':')[0]))

  pylab.figure(figsize=(10, 1.25*len(chrom_keys_sorted)))
  pylab.subplots_adjust(left=0.1, right=0.99, top=0.99, bottom=0.01, hspace=0.02)
  for n, key in enumerate(chrom_keys_sorted):
    # For each chromosome, plot a histogram of incorrectly aligned reads by position
    pylab.subplot(len(chrom_keys_sorted), 1, n+1)
    idx_wrong, = numpy.nonzero(~ra.correctly_aligned & (ra.correct_chrom == key))
    num_total = numpy.count_nonzero(ra.correct_chrom == key)
    seq_len = hdr[key]['len']
    if idx_wrong.size:
      pylab.hist(ra[idx_wrong].correct_pos, range=[0, seq_len], bins=100, normed=True, histtype='step')
    pylab.ylabel('{:s}\n({:2.2f}%)\n{:d}'.format(key, 100*float(idx_wrong.size)/num_total, num_total))
    pylab.setp(pylab.gca(), xlim=[0, seq_len], xticks=[], yticks=[])

  fig_fname = os.path.splitext(args['--file'])[0] + '_analyzed.pdf'
  pylab.savefig(fig_fname)
  print 'Saved to {:s}'.format(fig_fname)
