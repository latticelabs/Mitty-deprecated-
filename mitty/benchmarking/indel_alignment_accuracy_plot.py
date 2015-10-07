"""Plot alignment accuracy as a function of indel length. Accepts data in the .json format alindel puts out.
Plots one makes for paper, naturally, have to be customized, but this is a decent start.

  ----------------------------
  |A                         |
  |           Sim            |
  |         accuracy         |
  |                          |
  ----------------------------
  |B       read count        |
  ----------------------------
  |C        var count        |
  ----------------------------

"""
import itertools
from itertools import cycle
import cPickle
import os
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as ss
import click

mpl.rc('font', size=7)


def setup_figure():
  fig = plt.figure(figsize=(2*3.385, 2*3))  # two column figure for bio-informatics
  plt.subplots_adjust(left=0.15, bottom=0.1, right=0.98, top=0.93, wspace=0.05, hspace=0.01)
  gs = plt.GridSpec(3, 3,
                    width_ratios=[6, 0.9, 6],
                    height_ratios=[3.5, 1, 1])
  ax = {k: plt.subplot(g) for k, g in
        zip([''.join(e) for e in itertools.product(['A', 'B', 'C'], ['-DEL', '-REF/SNP', '-INS'])], gs)}
  return fig, ax


def decorate_axes(ax, axes_spec):
  """Decorate axes as needed"""
  # Set axes x-lims. We operate by columns here
  fxl = axes_spec['x-lim']
  x_lim = {'DEL': [fxl[0], -1], 'REF/SNP': [-1, 1], 'INS': [1, fxl[1]]}
  fxt = axes_spec['x-ticks']
  x_ticks = {'DEL': [t for t in fxt if t < 0], 'REF/SNP': [-0.5, 0.5], 'INS': [t for t in fxt if t > 0]}
  x_tick_labels = {'DEL': x_ticks['DEL'], 'REF/SNP': ['Ref', 'SNP'], 'INS': x_ticks['INS']}

  for k1 in ['DEL', 'REF/SNP', 'INS']:
    for k2 in ['A', 'B', 'C']:
      y_lim, y_ticks = axes_spec[k2]['y-lim'], axes_spec[k2]['y-ticks']
      y_tick_labels = axes_spec[k2].get('y-t-labels', y_ticks)
      ak = k2 + '-' + k1
      if k2 != 'A':
        ax[ak].set_yscale('log', nonposy='clip', subsy=[])
      plt.setp(ax[ak],
               xlim=x_lim[k1], xticks=x_ticks[k1], xticklabels=[] if k2 in ['A', 'B'] else x_tick_labels[k1],
               ylim=y_lim,
               yticks=[] if k1 != 'DEL' else y_ticks,
               yticklabels=[] if k1 != 'DEL' else y_tick_labels)
      if ak == 'C-REF/SNP':
        ax[ak].set_xticklabels(x_tick_labels[k1], rotation=90)

      ax[ak].xaxis.grid('true')  # The vertical guide lines

      # Ticks and axes
      for p in ['left', 'top', 'right']:
        if k1 == 'DEL' and p == 'left': continue
        ax[ak].spines[p].set_visible(False)

      ax[ak].minorticks_off()
      if k1 == 'DEL' or k2 == 'C':
        ax[ak].tick_params(direction='out', length=3, top='off', right='off', pad=0)
      else:
        ax[ak].tick_params(direction='in', left='off', top='off', right='off')

  ax['A-REF/SNP'].title.set(text=axes_spec['title'])
  ax['A-DEL'].legend(loc='lower left', fontsize=8)

  # X-labels
  ax['C-DEL'].set_xlabel('Deletions')
  ax['C-INS'].set_xlabel('Insertions')

  # Y-labels
  ax['A-DEL'].set_ylabel('% reads\ncorrectly aligned')
  ax['B-DEL'].set_ylabel('Read\ncount')
  ax['C-DEL'].set_ylabel('Variation\ncount')


def plot_indel_accuracy(cat_read_counts, pure_reference_read_counts, kernel_size=0, color='k', label=None, ax={}):
  """axs is a dict of axes A-DEL, A-REF/SNP and A-INS"""
  pc = cat_read_counts['correct'] / cat_read_counts['total'].astype(float) * 100
  if kernel_size > 2:
    pc = ss.medfilt(np.ma.masked_invalid(pc), kernel_size=kernel_size)

  idx = {
    'A-DEL': np.nonzero(cat_read_counts['x'] < 0)[0],
    'A-INS': np.nonzero(cat_read_counts['x'] > 0)[0],
  }
  for ak in ['A-DEL', 'A-INS']:
    ax[ak].plot(cat_read_counts['x'][idx[ak]], pc[idx[ak]], color=color, lw=1, label=label, alpha=0.71)

  # SNP
  idx = np.nonzero(cat_read_counts['x'] == 0)[0]
  x, y = 0.5, pc[idx]
  ax['A-REF/SNP'].plot(x, y, color + 'o', ms=2, mec=color, alpha=0.71)
  ax['A-REF/SNP'].plot(x, y, color + '_', lw=1.5, label=label, alpha=0.71)

  # ref reads
  ref_pc = pure_reference_read_counts[0] / float(pure_reference_read_counts[1]) * 100
  ax['A-REF/SNP'].plot(-0.5, ref_pc, color + 'o', ms=2, mec=color, alpha=0.71)
  ax['A-REF/SNP'].plot(-0.5, ref_pc, color + '_', lw=1.5, label=label, alpha=0.71)


def plot_read_count(cat_read_counts, pure_reference_read_counts, color='k', ax=[]):
  idx = {
    'B-DEL': np.nonzero(cat_read_counts['x'] < 0)[0],
    'B-INS': np.nonzero(cat_read_counts['x'] > 0)[0],
  }
  for ak in ['B-DEL', 'B-INS']:
    ax[ak].plot(cat_read_counts['x'][idx[ak]], cat_read_counts['total'][idx[ak]], color, lw=1)

  # SNP
  idx = np.nonzero(cat_read_counts['x'] == 0)[0]
  ax['B-REF/SNP'].plot(0.5, cat_read_counts['total'][idx], color + 'o', ms=2, mec=color)
  ax['B-REF/SNP'].plot(0.5, cat_read_counts['total'][idx], color + '_', lw=3)

  # ref reads
  cnt = pure_reference_read_counts[1]
  ax['B-REF/SNP'].plot(-0.5, cnt, color + 'o', ms=2, mec=color)
  ax['B-REF/SNP'].plot(-0.5, cnt, color + '_', lw=3)


def plot_var_count(indel_counts, color='k', ax=[]):
  idx = {
    'C-DEL': np.nonzero(indel_counts['x'] < 0)[0],
    'C-INS': np.nonzero(indel_counts['x'] > 0)[0],
  }
  for ak in ['C-DEL', 'C-INS']:
    ax[ak].plot(indel_counts['x'][idx[ak]], indel_counts['total'][idx[ak]], color, lw=1)

  # SNP
  idx = np.nonzero(indel_counts['x'] == 0)[0]
  ax['C-REF/SNP'].plot(0.5, indel_counts['total'][idx], color + 'o', ms=2, mec=color)
  ax['C-REF/SNP'].plot(0.5, indel_counts['total'][idx], color + '_', lw=3)


def load_indel_data(fname):
  data = json.load(open(fname, 'r'))
  data['indel_count'] = np.rec.fromrecords(data['indel_count'], dtype=[('x', 'int32'), ('total', 'uint32')])
  data['reads_within_feature'] = np.rec.fromrecords(data['reads_within_feature'], dtype=[('x', 'int32'), ('correct', 'uint32'), ('total', 'uint32')])
  data['templates_within_feature_but_read_outside'] = np.rec.fromrecords(data['templates_within_feature_but_read_outside'], dtype=[('x', 'int32'), ('correct', 'uint32'), ('total', 'uint32')])
  return data


@click.command()
@click.option('-f', multiple=True, type=click.Path(exists=True), help='Indel file(s) to plot')
@click.option('-o', default='indel_plot.pdf', type=click.Path(), help='Output file name')
@click.option('-l', multiple=True, help='Label to go with each file')
@click.option('--win', default=0, help='Size of median filter window to smooth plots')
@click.option('--indel-range', default=50, help='Maximum indel length to plot')
@click.option('--title', help='Title', default='Aligner accuracy')
def cli(f, o, l, win, indel_range, title):
  if len(f) == 0:
    return

  colors = ['k', 'r', 'b', 'g', 'y', 'c']
  #color_cycle = cycle(colors)
  plot_specs = [{'fname': fname, 'label': label or os.path.basename(fname), 'color_spec': color_spec}
                for fname, label, color_spec in itertools.izip_longest(f, l, colors) if fname is not None]

  fig, ax = setup_figure()
  max_read_cnt, max_indel_cnt = 0, 0
  for this_plot_spec in plot_specs:
    data = load_indel_data(this_plot_spec['fname'])
    max_read_cnt = max(max_read_cnt, data['reads_within_feature']['total'].max(), data['fully_outside_features'][1])
    max_indel_cnt = max(max_indel_cnt, data['indel_count']['total'].max())
    plot_indel_accuracy(data['reads_within_feature'], data['fully_outside_features'],
                        kernel_size=win, color=this_plot_spec['color_spec'], label=this_plot_spec['label'], ax=ax)
    plot_read_count(data['reads_within_feature'], data['fully_outside_features'], color='k', ax=ax)
    plot_var_count(data['indel_count'], color='k', ax=ax)

  axes_specs = {
    'title': title,
    'x-lim': [-indel_range, indel_range],
    'x-ticks': [-indel_range, 0, indel_range],
    'A': {'y-lim': [0, 105], 'y-ticks': [25, 50, 75, 95, 100]},
    'B': {'y-lim': [1, max_read_cnt * 2],
          'y-ticks': [1, max_read_cnt/100, max_read_cnt/10, max_read_cnt],
          'y-t-labels': [1, max_read_cnt/100, max_read_cnt/10, max_read_cnt]},
    'C': {'y-lim': [1, max_indel_cnt * 2],
          'y-ticks': [1, max_indel_cnt/100, max_indel_cnt/10, max_indel_cnt],
          'y-t-labels': [1, max_indel_cnt/100, max_indel_cnt/10, max_indel_cnt]}
  }
  decorate_axes(ax, axes_specs)

  plt.savefig(o)

if __name__ == '__main__':
  cli()