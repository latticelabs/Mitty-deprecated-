"""Figure 1 of aligner paper

Panel 1:
x-axis No variant
y-axis Alignment accuracy

Panel 2:
x-axis length of variant
y-axis alignment accuracy

Algorithm:

1. Load the variant database
2. Load the categorized reads
3. For each variant, around a window, compute the total number of reads and the correct reads
4. Split the variants into bins and plot the data

"""
#from itertools import izip
import cPickle

import click
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import scipy.signal as ss


from scipy.stats import beta


# inspired by https://gist.github.com/paulgb/6627336
def pc_binom_interval(correct, total, confint=0.95):
    quantile = (1 - confint) / 2.
    lower = beta.ppf(quantile, correct + 1, total - correct + 1)
    upper = beta.ppf(1 - quantile, correct + 1, total - correct + 1)
    return lower, upper


def plot_indel_accuracy(cat_counts, kernel_size=0, color='k'):
  cat_counts = cat_counts[:-1]
  pc = cat_counts['correct'] / cat_counts['total'].astype(float) * 100
  #pc_ref = ref_read_counts[0, 0] / float(ref_read_counts[0, 1]) * 100
  correct, total = cat_counts['correct'], cat_counts['total']
  if kernel_size > 2:
    pc = ss.medfilt(np.ma.masked_invalid(pc), kernel_size=kernel_size)
    correct = pc * total
  l, u = pc_binom_interval(correct, total)
  if kernel_size > 2:
    l = ss.medfilt(np.ma.masked_invalid(l), kernel_size=kernel_size)
    u = ss.medfilt(np.ma.masked_invalid(u), kernel_size=kernel_size)

  plt.fill_between(cat_counts['x'], l, u, color='', facecolor=color, alpha=.1, interpolate=True)
  plt.plot(cat_counts['x'], pc, color, lw=2)


def adjust_spines(ax, spines):
  for loc, spine in ax.spines.items():
    if loc in spines:
      spine.set_position(('outward', 10))  # outward by 10 points
      spine.set_smart_bounds(True)
    else:
      spine.set_color('none')  # don't draw spine

  # turn off ticks where there is no spine
  if 'left' in spines:
    ax.yaxis.set_ticks_position('left')
  else:
    # no yaxis ticks
    ax.yaxis.set_ticks([])

  if 'bottom' in spines:
    ax.xaxis.set_ticks_position('bottom')
  else:
    # no xaxis ticks
    ax.xaxis.set_ticks([])


@click.command()
@click.option('-f', multiple=True, type=click.Path(exists=True), help='Indel file(s) to plot')
@click.option('-o', default='indel_plot.pdf', type=click.Path(), help='Output file name')
@click.option('-l', multiple=True, help='Label to go with each file')
@click.option('--win', default=0, help='Size of median filter window to smooth plots')
@click.option('--indel-range', default=500, help='Maximum indel length to plot')
def cli(f, o, l, win, indel_range):
  if len(f) == 0:
    return

  if len(l) < len(f):
    l = [str(n) for n in range(len(f))]

  fig = plt.figure()

  ylim = [0, 105]
  yticks = [25, 50, 75, 95, 96, 97, 98, 99, 100]

  plt.subplots_adjust(bottom=0.22, top=0.95, left=0.1, right=0.98, wspace=0.05)
  ax2 = plt.subplot2grid((1, 4), (0, 1), colspan=3)
  colors = ['k', 'r', 'b', 'c', 'm', 'y', 'g']

  pc_ref_var = []
  pc_ref_ref = []
  for fn, col in zip(f, colors[:len(f)]):
    d = cPickle.load(open(fn))
    plot_indel_accuracy(d['cc'], color=col, kernel_size=win)
    pc_ref_var.append(d['rr'][0, 0] / float(d['rr'][0, 1]) * 100)
    pc_ref_ref.append(d['rr'][1, 0] / float(d['rr'][1, 1]) * 100)

  plt.setp(ax2, yticks=yticks, yticklabels=[], ylim=ylim, xlim=[-indel_range, indel_range])
  ax2.spines['right'].set_visible(False)
  ax2.spines['top'].set_visible(False)
  plt.xlabel('Indel length')
  plt.title('Variant regions')

  plt.plot([0, 0], [15, 105], 'k--')
  plt.text(0, 5, 'SNP', horizontalalignment='center')
  plt.text(300, 5, 'Insertions', horizontalalignment='center')
  plt.text(-300, 5, 'Deletions', horizontalalignment='center')

  ax1 = plt.subplot2grid((1, 4), (0, 0), sharey=ax2)
  x = np.arange(len(pc_ref_ref))
  width = 0.3
  ax1.bar(x - width, pc_ref_var, width, color='y')
  ax1.bar(x, pc_ref_ref, width, color='c')
  plt.setp(ax1, yticks=yticks, ylim=ylim, xlim=[-1, len(pc_ref_ref)])
  plt.xticks(x, l, rotation='vertical')
  ax1.spines['right'].set_visible(False)
  ax1.spines['top'].set_visible(False)
  plt.ylabel('% correct')
  plt.title('Null regions')

  plt.savefig(o)

if __name__ == '__main__':
  cli()