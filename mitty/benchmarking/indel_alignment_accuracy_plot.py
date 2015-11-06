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


def setup_figure(diff=False):
  """Set diff to True if you want an additional panel showing pair-wise differences in accuracy"""
  fig = plt.figure(figsize=(2*3.385, 2*3))  # two column figure for bio-informatics
  plt.subplots_adjust(left=0.15, bottom=0.1, right=0.98, top=0.93, wspace=0.05, hspace=0.01)
  gs = plt.GridSpec(4 if diff else 3, 3,
                    width_ratios=[6, 0.9, 6],
                    height_ratios=[3.5, 2, 1, 1] if diff else [3.5, 1, 1])
  ax = {k: plt.subplot(g) for k, g in
        zip([''.join(e) for e in itertools.product(['A', 'Ad', 'B', 'C'] if diff else ['A', 'B', 'C'], ['-DEL', '-REF/SNP', '-INS'])], gs)}
  return fig, ax


def decorate_axes(ax, axes_spec):
  """Decorate axes as needed"""
  # Set axes x-lims. We operate by columns here
  fxl = axes_spec['x-lim']
  x_lim = {'DEL': [fxl[0], -1], 'REF/SNP': [-1, 1], 'INS': [1, fxl[1]]}
  fxt = axes_spec['x-ticks']
  x_ticks = {'DEL': [t for t in fxt if t < 0], 'REF/SNP': [-0.5, 0.5], 'INS': [t for t in fxt if t > 0]}
  x_tick_labels = {'DEL': [abs(x) for x in x_ticks['DEL']], 'REF/SNP': ['Ref', 'SNP'], 'INS': x_ticks['INS']}

  for k1 in ['DEL', 'REF/SNP', 'INS']:
    for k2 in ['A', 'Ad', 'B', 'C']:
      if k2 not in axes_spec: continue
      y_lim, y_ticks = axes_spec[k2]['y-lim'], axes_spec[k2]['y-ticks']
      y_tick_labels = axes_spec[k2].get('y-t-labels', y_ticks)
      ak = k2 + '-' + k1
      if k2 not in ['A', 'Ad']:
        ax[ak].set_yscale('log', nonposy='clip', subsy=[])
      plt.setp(ax[ak],
               xlim=x_lim[k1], xticks=x_ticks[k1], xticklabels=[] if k2 != 'C' else x_tick_labels[k1],
               ylim=y_lim,
               yticks=y_ticks,
               yticklabels=[] if k1 != 'DEL' else y_tick_labels)
      if ak == 'C-REF/SNP':
        ax[ak].set_xticklabels(x_tick_labels[k1], rotation=90)

      ax[ak].xaxis.grid('true')  # The vertical guide lines

      if k2 in ['A', 'Ad']:
        ax[ak].yaxis.grid('true')  # The horizontal guide lines

      # Ticks and axes
      for p in ['left', 'top', 'right']:
        if k1 == 'DEL' and p == 'left': continue
        ax[ak].spines[p].set_visible(False)

      ax[ak].minorticks_off()
      ax[ak].tick_params(bottom='off', left='off', top='off', right='off')
      if k2 == 'C':
        ax[ak].tick_params(direction='out', length=3, bottom='on', pad=3)
      if k1 == 'DEL':
        ax[ak].tick_params(direction='out', length=3, left='on', pad=2)

  ax['A-REF/SNP'].title.set(text=axes_spec['title'])
  ax['A-DEL'].legend(loc='lower left', fontsize=8)
  if 'Ad-DEL' in ax: ax['Ad-DEL'].legend(loc='lower left', fontsize=8)

  # X-labels
  ax['C-DEL'].set_xlabel('Deletions (bp)')
  ax['C-INS'].set_xlabel('Insertions (bp)')

  # Y-labels
  ax['A-DEL'].set_ylabel('% reads\ncorrectly aligned')
  if 'Ad-DEL' in ax: ax['Ad-DEL'].set_ylabel('Difference')
  ax['B-DEL'].set_ylabel('Read\ncount')
  ax['C-DEL'].set_ylabel('Variation\ncount')


#def plot_indel_accuracy(cat_read_counts, pure_reference_read_counts, kernel_size=0, color='k', label=None, ax={}):
def plot_indel_accuracy(indel_size, pc_raw, pc, ref_pc, color='k', label=None, ax={}, plot_row='A'):
  """axs is a dict of axes A-DEL, A-REF/SNP and A-INS"""
  #pc_raw = cat_read_counts['correct'] / cat_read_counts['total'].astype(float) * 100
  #pc = ss.medfilt(np.ma.masked_invalid(pc_raw), kernel_size=kernel_size) if kernel_size > 2 else pc_raw

  idx = {
    plot_row + '-DEL': np.nonzero(indel_size < 0)[0],
    plot_row + '-INS': np.nonzero(indel_size > 0)[0],
  }
  for ak in [plot_row + '-DEL', plot_row + '-INS']:
    ax[ak].plot(indel_size[idx[ak]], pc[idx[ak]], color=color, lw=1, label=label, alpha=0.71)

  # SNP
  idx = np.nonzero(indel_size == 0)[0]
  x, y = 0.5, pc_raw[idx]
  ax[plot_row + '-REF/SNP'].plot(x, y, color + 'o', ms=2, mec=color, alpha=0.71)
  ax[plot_row + '-REF/SNP'].plot(x, y, color + '_', lw=1.5, label=label, alpha=0.71)

  # ref reads
  #ref_pc = pure_reference_read_counts[0] / float(pure_reference_read_counts[1]) * 100
  ax[plot_row + '-REF/SNP'].plot(-0.5, ref_pc, color + 'o', ms=2, mec=color, alpha=0.71)
  ax[plot_row + '-REF/SNP'].plot(-0.5, ref_pc, color + '_', lw=1.5, label=label, alpha=0.71)


def plot_read_count(indel_size, indel_reads, ref_reads, color='k', ax=[]):
  idx = {
    'B-DEL': np.nonzero(indel_size < 0)[0],
    'B-INS': np.nonzero(indel_size > 0)[0],
  }
  for ak in ['B-DEL', 'B-INS']:
    ax[ak].plot(indel_size[idx[ak]], indel_reads[idx[ak]], color, lw=1)

  # SNP
  idx = np.nonzero(indel_size == 0)[0]
  ax['B-REF/SNP'].plot(0.5, indel_reads[idx], color + 'o', ms=2, mec=color)
  ax['B-REF/SNP'].plot(0.5, indel_reads[idx], color + '_', lw=3)

  # ref reads
  ax['B-REF/SNP'].plot(-0.5, ref_reads, color + 'o', ms=2, mec=color)
  ax['B-REF/SNP'].plot(-0.5, ref_reads, color + '_', lw=3)


def plot_var_count(indel_size, indel_counts, color='k', ax=[]):
  idx = {
    'C-DEL': np.nonzero(indel_size < 0)[0],
    'C-INS': np.nonzero(indel_size > 0)[0],
  }
  for ak in ['C-DEL', 'C-INS']:
    ax[ak].plot(indel_size[idx[ak]], indel_counts[idx[ak]], color, lw=1)

  # SNP
  idx = np.nonzero(indel_size == 0)[0]
  ax['C-REF/SNP'].plot(0.5, indel_counts[idx], color + 'o', ms=2, mec=color)
  ax['C-REF/SNP'].plot(0.5, indel_counts[idx], color + '_', lw=3)


def load_indel_data(fname, label):
  data = json.load(open(fname, 'r'))
  data['indel_count'] = np.rec.fromrecords(data['indel_count'], dtype=[('x', 'int32'), ('total', 'uint32')])
  data['reads_within_feature'] = np.rec.fromrecords(data['reads_within_feature'], dtype=[('x', 'int32'), ('correct', 'uint32'), ('total', 'uint32')])
  data['templates_within_feature_but_read_outside'] = np.rec.fromrecords(data['templates_within_feature_but_read_outside'], dtype=[('x', 'int32'), ('correct', 'uint32'), ('total', 'uint32')])
  data['label'] = label
  return data


def process_indel_data(data, win):
  _pc_raw = data['reads_within_feature']['correct'] / data['reads_within_feature']['total'].astype(float) * 100
  pc_raw = np.ma.masked_invalid(_pc_raw)
  pc = ss.medfilt(pc_raw, kernel_size=win) if win > 2 else pc_raw
  ref_pc = data['fully_outside_features'][0] / float(data['fully_outside_features'][1]) * 100
  indel_read_total = data['reads_within_feature']['total']
  ref_read_total = data['fully_outside_features'][1]
  return {
    'label': data['label'],
    'indel_count': data['indel_count']['total'],
    'indel_size': data['indel_count']['x'],
    'pc_raw': pc_raw,
    'pc': pc,
    'ref_pc': ref_pc,
    'indel_read_total': indel_read_total,
    'ref_read_total': ref_read_total
  }


def process_pairwise_data(all_data, i, j):
  return {
    'label': all_data[i]['label'] + '-' + all_data[j]['label'],
    'indel_count': all_data[i]['indel_count'],
    'indel_size': all_data[i]['indel_size'],  # We assume this is identical for all runs!
    'pc_raw': all_data[i]['pc_raw'] - all_data[j]['pc_raw'],
    'pc': all_data[i]['pc'] - all_data[j]['pc'],
    'ref_pc': all_data[i]['ref_pc'] - all_data[j]['ref_pc']
  }


def build_axes_specs(all_data, pairwise_data, indel_range, title):
  dta = all_data[0]
  max_read_cnt, min_read_cnt = max(2, dta['ref_read_total'], dta['indel_read_total'].max()), min(1, dta['ref_read_total'], dta['indel_read_total'].min())
  max_indel_cnt, min_indel_cnt = max(2, dta['indel_count'].max()), min(1, dta['indel_count'].min())
  if indel_range is None:
    indel_range = max([abs(n) for n in dta['indel_size']])
  axes_specs = {
    'title': title,
    'x-lim': [-indel_range, indel_range],
    'x-ticks': [-indel_range, -1, 0, 1, indel_range],
    'A': {'y-lim': [0, 105], 'y-ticks': [25, 50, 75, 95, 100]},
    'B': {'y-lim': [min_read_cnt / 2, max_read_cnt * 2],
          'y-ticks': [min_read_cnt, max_read_cnt],
          'y-t-labels': [min_read_cnt, max_read_cnt]},
    'C': {'y-lim': [min_indel_cnt / 2, max_indel_cnt * 2],
          'y-ticks': [min_indel_cnt, max_indel_cnt],
          'y-t-labels': [min_indel_cnt, max_indel_cnt]}
  }
  if len(pairwise_data) > 0:
    max_delta_pc = max(0, max([max(np.ma.max(pwd['pc']), pwd['pc_raw'][np.nonzero(pwd['indel_size'] == 0)[0]]) for pwd in pairwise_data]))
    min_delta_pc = min(0, min([min(np.ma.min(pwd['pc']), pwd['pc_raw'][np.nonzero(pwd['indel_size'] == 0)[0]]) for pwd in pairwise_data]))
    axes_specs['Ad'] = {'y-lim': [min_delta_pc * 1.05, max_delta_pc * 1.05],
                        'y-ticks': [min_delta_pc, 0, max_delta_pc],
                        'y-t-labels': [round(min_delta_pc, 3), 0, round(max_delta_pc, 3)]}
  return axes_specs


@click.command()
@click.option('-f', multiple=True, type=click.Path(exists=True), help='Indel file(s) to plot')
@click.option('-o', default='indel_plot.pdf', type=click.Path(), help='Output file name')
@click.option('-l', multiple=True, help='Label to go with each file')
@click.option('--win', default=0, help='Size of median filter window to smooth plots')
@click.option('--indel-range', type=int, default=None, help='Maximum indel length to plot')
@click.option('--title', help='Title', default='Aligner accuracy')
def cli(f, o, l, win, indel_range, title):
  if len(f) == 0:
    return

  colors = ['k', 'r', 'b', 'g', 'y', 'c']
  #color_cycle = cycle(colors)
  plot_specs = [{'fname': fname, 'label': label or os.path.basename(fname), 'color_spec': color_spec}
                for fname, label, color_spec in itertools.izip_longest(f, l, colors) if fname is not None]
  all_data = [process_indel_data(load_indel_data(this_plot_spec['fname'], this_plot_spec['label']), win) for this_plot_spec in plot_specs]
  pairwise_data = [process_pairwise_data(all_data, i, j) for i in range(len(all_data) - 1) for j in range(i + 1, len(all_data))]

  axes_specs = build_axes_specs(all_data, pairwise_data, indel_range, title)

  fig, ax = setup_figure(diff='Ad' in axes_specs)

  for this_plot_spec, data in itertools.izip(plot_specs, all_data):
    plot_indel_accuracy(data['indel_size'], data['pc_raw'], data['pc'], data['ref_pc'],
                        color=this_plot_spec['color_spec'], label=data['label'], ax=ax,
                        plot_row='A')
  if len(all_data):
    dta = all_data[0]
    plot_read_count(dta['indel_size'], dta['indel_read_total'], dta['ref_read_total'], color='k', ax=ax)
    plot_var_count(dta['indel_size'], dta['indel_count'], color='k', ax=ax)

  for clr_ctr, pwd in enumerate(pairwise_data):
    plot_indel_accuracy(pwd['indel_size'], pwd['pc_raw'], pwd['pc'], pwd['ref_pc'],
                        color=colors[clr_ctr], label=pwd['label'], ax=ax, plot_row='Ad')

  decorate_axes(ax, axes_specs)
  plt.savefig(o)

if __name__ == '__main__':
  cli()