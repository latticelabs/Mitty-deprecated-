"""Take in the output of cheata split and produce a Python pickle file of alignment diagnostics and a summary figure
for the sequence. This code also requires the .smalla file of the original sequence to plot the complexity measure
across the sequence.

The code provides the following data (and summarizes it on one page)

1. List of correct and actual positions of incorrectly placed reads
2. Complexity score of the sequence of each read
3. Mapping quality of each read as reported by the aligner

Usage:
diagnose_alignment.py --correctbam=CBAM --wrongbam=WBAM  --smalla=SMALLA  --outpkl=OUTPKL  --outfig=OUTFIG
diagnose_alignment.py test [-v]

Options:
  --correctbam=CBAM    -  .bam file (with index) of reads correctly placed by the aligner
  --wrongbam=WBAM      -  .bam file (with index) of reads incorrectly placed by the aligner
  --smalla=SMALLA      -  .smalla file of the orginal reference sequence (needed for sequence entropy plot)
  --outpkl=OUTPKL      -  .pkl file with alignment diagnostics
  --outfig=OUTFIG      -  .pdf/.png/.jpg file of the summary figure of alignment diagnostics
"""
import matplotlib
matplotlib.use('Agg')
import math
import pysam
import pylab
import mmap
import docopt
import os
import cPickle

# TODO: figure out if a function call in a list comprehension is cheaper than a function call in a loop

def analyse_bam(bam):
  """Given an input bam file, for each read compute the following:
    1. List of correct and actual position
    2. Complexity score of the sequence
    3. Mapping quality as reported by the aligner
  """
  bam_file = pysam.Samfile(bam, 'rb')
  hdr = bam_file.header
  aligned_pos = []
  correct_pos = []
  mapping_score = []
  sequence_entropy = []
  for read in bam_file:
    cheat_answer = read.qname.split(':')
    this_correct_pos = int(cheat_answer[1])
    if read.flag & 0x01:  # Paired reads
      if read.flag & 0x80:  # Second end (mate)
        this_correct_pos = int(cheat_answer[3])

    computed_pos = read.pos + 1  # PySam uses 0 indexing ...
    aligned_pos.append(computed_pos)
    correct_pos.append(this_correct_pos)
    mapping_score.append(read.mapq)
    sequence_entropy.append(shannon_entropy(read.seq))
  bam_file.close()
  return {
    'bam name': bam,
    'sequence name': hdr['SQ'][0]['SN'],
    'sequence len': hdr['SQ'][0]['LN'],
    'read count': len(aligned_pos),
    'aligned pos': aligned_pos,
    'correct pos': correct_pos,
    'mapping score': mapping_score,
    'sequence entropy': sequence_entropy
  }


def shannon_entropy(sequence, alphabet=['A', 'C', 'T', 'G']):
  """Given a sequence compute the shannon entropy given the alphabet. Do the computation in blocks if asked for
  >>> shannon_entropy('AACCTTGG')
  2.0
  """
  seq_len = float(len(sequence))
  p_obs_alphabet = [sequence.count(k) / seq_len for k in alphabet]
  return - sum([p_obs * math.log(p_obs, 2) for p_obs in p_obs_alphabet if p_obs > 0])


def shannon_entropy_blocked(sequence, block_size=10000, block_overlap=5000, alphabet=['A', 'C', 'T', 'G']):
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
  fig = pylab.figure(figsize=(9, 15))
  fig.subplots_adjust(top=0.97, bottom=0.1, left=0.09, right=0.99)
  H = {
    'scatter': pylab.subplot2grid((8, 2), (0, 0), colspan=2, rowspan=5, xticks=[]),
    'alignment error': pylab.subplot2grid((8, 2), (5, 0), colspan=2, xticks=[], yticks=[]),
    'entropy': pylab.subplot2grid((8, 2), (6, 0), colspan=2, xticks=[], yticks=[]),
    'map score hist': pylab.subplot2grid((8, 2), (7, 0)),
    'local entropy hist': pylab.subplot2grid((8, 2), (7, 1))
  }
  return H


def alignment_error_scatter_plot(bad):
  x = pylab.array(bad['correct pos'])
  y = pylab.array(bad['aligned pos'])
  ms = pylab.array(bad['mapping score'])
  idx_lo = pylab.find(ms < 40)
  pylab.plot(x[idx_lo], y[idx_lo], 'k.')
  idx_hi = pylab.find(ms > 40)
  pylab.plot(x[idx_hi], y[idx_hi], 'yo')
  pylab.plot([0, bad['sequence len']], [0, bad['sequence len']], 'k:')
  pylab.axis('scaled')
  pylab.setp(pylab.gca(), xlim=[0, bad['sequence len']], ylim=[0, bad['sequence len']], xlabel='Correct pos', ylabel='Aligned pos')


def alignment_error_density_plot(bad):
  if len(bad['correct pos']):
    pylab.hist(bad['correct pos'], range=[0, bad['sequence len']], bins=1000, normed=True, histtype='step', lw=3)
  pylab.setp(pylab.gca(), xlim=[0, bad['sequence len']])
  pylab.ylabel('Error density')


def entropy_plot(seq):
  entropy, index = shannon_entropy_blocked(seq, block_size=1000, block_overlap=0, alphabet=['A', 'C', 'T', 'G'])
  pylab.plot(index, entropy)
  pylab.setp(pylab.gca(), xlim=[0, len(seq)], ylim=[0, 2.1], ylabel='Bits')
  #pylab.gca().set_yscale("log", nonposy='clip')


def alignment_score_histogram(good, bad):
  if len(bad['mapping score']):
    pylab.hist(bad['mapping score'], bins=21, range=[0,70], color='r', histtype='step', lw=3)
  if len(good['mapping score']):
    pylab.hist(good['mapping score'], bins=21, range=[0, 70], color='y', histtype='step', lw=2)
  pylab.gca().set_yscale("log", nonposy='clip')
  pylab.xlabel('Mapping score')
  pylab.ylabel('Read count')


def local_entropy_histogram(good, bad):
  if len(bad['sequence entropy']):
    pylab.hist(bad['sequence entropy'], range=[0, 2], bins=1000, normed=True, histtype='step', lw=3, color='r')
  if len(good['sequence entropy']):
    pylab.hist(good['sequence entropy'], range=[0, 2], bins=1000, normed=True, histtype='step', lw=1, color='g')
  pylab.gca().set_yscale("log", nonposy='clip')
  pylab.xlabel('Local entropy')
  pylab.ylabel('Read density')


def main(goodbam, badbam, smalla_fname, pkl_fname, out_fig_fname):
  f_seq = open(smalla_fname, 'rb')
  seq = mmap.mmap(f_seq.fileno(), 0, access=mmap.ACCESS_READ)

  good = analyse_bam(goodbam)
  bad = analyse_bam(badbam)

  H = setup_axes()

  pylab.axes(H['scatter'])
  alignment_error_scatter_plot(bad)

  pylab.axes(H['alignment error'])
  alignment_error_density_plot(bad)

  pylab.axes(H['entropy'])
  entropy_plot(seq)

  pylab.axes(H['map score hist'])
  alignment_score_histogram(good, bad)

  pylab.axes(H['local entropy hist'])
  local_entropy_histogram(good, bad)

  title = '{:s}: {:0.2f}% of {:d} reads misaligned'.format(os.path.basename(goodbam).split('_')[0],
                                                 100 * bad['read count']/float(bad['read count'] + good['read count']),
                                                 bad['read count'] + good['read count'])

  pylab.suptitle(title)
  pylab.savefig(out_fig_fname)

  with open(pkl_fname,'w') as f_pkl:
    cPickle.dump({'good': good, 'bad': bad}, f_pkl, protocol=cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__)
  if args['test']:
    import doctest
    doctest.testmod()
    exit(0)

  main(args['--correctbam'], args['--wrongbam'], args['--smalla'], args['--outpkl'], args['--outfig'])
