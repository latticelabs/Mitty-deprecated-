"""Methods and routines useful for the mappability project.

Usage:
mappability.py <pklfname> <smallafname> <outfigfname>
mappability.py test

Options:

<pklfname>        -   name of the .pkl file from cheata check
<smallafname>     -   name of the smalla file used as reference
<outfigfname>     -   output figure file name
"""
import matplotlib
matplotlib.use('Agg')
import math
import cPickle
import pylab
import mmap
import docopt
import os

def shannon_entropy(sequence, block_size=10000, block_overlap=5000, alphabet=['A', 'C', 'T', 'G']):
  """Given a sequence compute the shannon entropy given the alphabet. Do the computation in blocks if asked for

  >>> shannon_entropy('AACCTTGG')
  ([2.0], [4.0])
  >>> shannon_entropy('AACCTTGG'*8, block_size=8, block_overlap=1)
  ([2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
   [4.0, 11.0, 18.0, 25.0, 32.0, 39.0, 46.0, 53.0, 60.0])
  """
  seq_end = len(sequence)
  blk_advance = block_size - block_overlap
  assert blk_advance > 0
  blk_start = 0
  entropy = []
  index = []
  while blk_start < seq_end - 1:
    blk_end = min(blk_start + block_size, seq_end)
    sub_sequence = sequence[blk_start:blk_end]
    sub_sequence_len = float(blk_end - blk_start)
    p_obs_alphabet = [sub_sequence.count(k) / sub_sequence_len for k in alphabet]
    entropy.append(- sum([p_obs * math.log(p_obs, 2) for p_obs in p_obs_alphabet if p_obs > 0]))
    index.append((blk_start + blk_end)/2.0)
    blk_start += blk_advance
  return entropy, index


def setup_axes():
  from matplotlib.patches import Rectangle as pRect
  fig = pylab.figure(figsize=(9.8,14))
  H = {
    'scatter': pylab.axes([0.07, .29, .92, .67], xticks=[]),
    'alignment error': pylab.axes([0.07, .16, .92, .11], xticks=[]),
    'entropy': pylab.axes([0.07, .04, .92, .11]),
    'map score hist': pylab.axes([.15, .75, .3, .12], axisbg=None)
  }
  H['map score hist'].patch.set_alpha(0.5)
  rect = pRect((0.1, 0.73), 0.32, 0.13, facecolor='white', edgecolor='none', transform=fig.transFigure, zorder=-1)
  fig.patches.append(rect)

  return H


def alignment_error_scatter_plot(st):
  x = st['correct_pos']
  y = st['aligned_pos']
  pylab.plot(x,y,'k.')
  pylab.plot([0, st['len']], [0, st['len']], 'k:')
  pylab.axis('scaled')
  pylab.setp(pylab.gca(), xlim=[0, st['len']], ylim=[0, st['len']], xlabel='Correct pos', ylabel='Aligned pos')


def alignment_score_histogram(st):
  pylab.hist(st['bad_alignment_map_score'], bins=21, range=[0,40], color='r', histtype='step', lw=5)
  pylab.hist(st['good_alignment_map_score'], bins=21, range=[0,40], color='y', histtype='step', lw=5)
  pylab.gca().set_yscale("log", nonposy='clip')
  pylab.xlabel('Mapping score')
  pylab.ylabel('Read count')


def entropy_plot(seq):
  entropy, index = shannon_entropy(seq, block_size=10000, block_overlap=0, alphabet=['A', 'C', 'T', 'G'])
  pylab.plot(index, entropy)
  pylab.setp(pylab.gca(), xlim=[0, len(seq)], ylim=[0, 2], ylabel='Bits')


def alignment_error_density_plot(st):
  pylab.hist(st['correct_pos'], range=[0, st['len']], bins=1000, normed=True, histtype='step')
  pylab.setp(pylab.gca(), xlim=[0, st['len']])
  pylab.ylabel('Error density')


def main(pkl_fname, smalla_fname):
  with open(pkl_fname,'r') as f_pkl, open(smalla_fname, 'rb') as f_seq:
    st = cPickle.load(f_pkl)
    seq = mmap.mmap(f_seq.fileno(), 0, access=mmap.ACCESS_READ)

    H = setup_axes()

    pylab.axes(H['scatter'])
    alignment_error_scatter_plot(st)

    pylab.axes(H['alignment error'])
    alignment_error_density_plot(st)

    pylab.axes(H['entropy'])
    entropy_plot(seq)

    pylab.axes(H['map score hist'])
    alignment_score_histogram(st)

    misaligned_reads = len(st['bad_alignment_map_score'])
    total_reads = misaligned_reads + len(st['good_alignment_map_score'])
    title='{:s}: {:0.2f}% reads misaligned'.format(os.path.basename(args['<pklfname>']).split('_')[0],
                                                   100 * misaligned_reads / float(total_reads))
    pylab.suptitle(title)

if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__)
  if args['test']:
    import doctest
    doctest.testmod()
    exit(0)

  main(args['<pklfname>'], args['<smallafname>'])
  pylab.savefig(args['<outfigfname>'])