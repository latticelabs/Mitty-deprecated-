"""Plot misalignments as circular plots or 2D histograms

Commandline::

  Usage:
    alignment (circle|matrix) <bam_file> [--down-sample=S] [pdf]

  Options:
    circle           Plot a circle plot
    matrix           Plot a 2D histogram
    <bam_file>       Name of original bam file
    --down-sample=S  Factor by which to down-sample
    pdf              Should the output be pdf? (png otherwise)

"""
import csv
import docopt
import json
import os

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
import numpy as np


def draw_circle_plot(chrom_lens, mis):
  ax1 = plt.subplot(111, polar=True)
  plot_genome_as_a_circle(chrom_lens, radius=1, chrom_gap=0.05, chrom_thick=5)
  plot_mis_alignments_on_a_circle(chrom_lens, mis, radius=1, chrom_gap=0.05, lw=0.5)


def plot_genome_as_a_circle(chrom_lens, radius=1, chrom_gap=0.001, chrom_thick=0.1):
  """Plot the chromosomes on a circle. Will plot on the currently active axes. Needs axes to be polar

  :param (list) chrom_lens: list of chromosome lengths
  :param (float) radius: radius of the plot
  :param (float) chrom_gap: gap between chromosomes in radians
  :param (float) chrom_thick: thickness of chromosome
  """
  total_len = sum(chrom_lens)
  radians_per_base = (2.0 * np.pi - len(chrom_lens) * chrom_gap) / total_len  # With allowance for chrom gaps

  xticks = []
  xticklabels = []
  delta_radian = 0.01
  start_radian = 0
  for ch_no, l in enumerate(chrom_lens):
    end_radian = start_radian + l * radians_per_base
    theta = np.arange(start_radian, end_radian, delta_radian)
    plt.plot(theta, [radius] * theta.size, lw=chrom_thick)
    xticks.append((start_radian + end_radian)/2)
    xticklabels.append(str(ch_no + 1))
    start_radian = end_radian + chrom_gap

  ax = plt.gca()
  plt.setp(ax.get_yticklabels(), visible=False)
  ax.grid(False)
  plt.setp(ax, xticks=xticks, xticklabels=xticklabels)


def plot_mis_alignments_on_a_circle(chrom_lens, misalignments, radius=1, chrom_gap=0.001, lw=2):
  """.

  :param (list) chrom_lens: list of chromosome lengths
  :param (list) misalignments: list of tuples (correct_chrom, correct_pos, aligned_chrom, aligned_pos)
  :param (float) radius: radius of the plot
  :param (float) lw: thickness of drawn line
  """
  total_len = sum(chrom_lens)
  radians_per_base = (2.0 * np.pi - len(chrom_lens) * chrom_gap) / total_len  # With allowance for chrom gaps

  # http://matplotlib.org/users/path_tutorial.html
  codes = [
    Path.MOVETO,
    Path.CURVE4,
    Path.CURVE4,
    Path.CURVE4,
  ]
  ax = plt.gca()

  for m in misalignments:
    t0 = (sum(chrom_lens[:m[0]-1]) + m[1]) * radians_per_base + (m[0]-1) * chrom_gap
    t1 = (sum(chrom_lens[:m[2]-1]) + m[3]) * radians_per_base + (m[2]-1) * chrom_gap
    verts = [
      (t0, radius),       # P0
      (t0, .2 * radius),  # P1
      (t1, .2 * radius),  # P2
      (t1, radius),       # P3
    ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', lw=lw)
    ax.add_patch(patch)
    plt.plot(t0, .99 * radius, 'x')


def draw_matrix_plot(chrom_lens, mis):
  flattened_coordinates, chrom_offsets = flatten_coordinates(chrom_lens, mis)
  plot_mis_alignment_matrix(flattened_coordinates, chrom_offsets, histogram=False)


def flatten_coordinates(chrom_lens, mis):
  # First work out the chromosome offsets
  chrom_offsets = reduce(lambda x, y: x + [x[-1] + y], chrom_lens, [0])
  flattened_coordinates = np.array([[chrom_offsets[r[0] - 1] + r[1], chrom_offsets[r[2] - 1] + r[3]] for r in mis])
  return flattened_coordinates, chrom_offsets


def plot_mis_alignment_matrix(data, chrom_offsets, histogram=False):
  """
  :param (ndarray) data: n x 2 array of read positions (unrolled to whole genome) 1st column = correct 2nd column computed
  """
  if histogram:
    h, xe, ye = np.histogram2d(data[:, 1], data[:, 0], bins=100, range=[[0, chrom_offsets[-1]], [0, chrom_offsets[-1]]])
    #plt.imshow(h, origin="lower", extent=(0, chrom_offsets[-1], 0, chrom_offsets[-1]), cmap=plt.cm.gray_r, norm=LogNorm(), interpolation=None)
    plt.matshow(h, origin="lower", extent=(0, chrom_offsets[-1], 0, chrom_offsets[-1]), cmap=plt.cm.gray_r, norm=LogNorm(), interpolation=None)
  else:
    plt.scatter(data[:,0], data[:,1], s=1, marker='.')
    ax = plt.gca()
    plt.setp(ax, xlim=[0, chrom_offsets[-1]], ylim=[0, chrom_offsets[-1]], aspect=1)

  ax = plt.gca()
  plt.setp(ax.get_xticklabels(), visible=False)
  plt.setp(ax.get_yticklabels(), visible=False)
  plt.setp(ax, xticks=chrom_offsets, yticks=chrom_offsets, xlabel='Correct position', ylabel='Aligned position')


def cli():
  """Command line script entry point."""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__)

  prefix = os.path.splitext(args['<bam_file>'])[0]
  summary_fname = prefix + '_summary.json'
  misalignments_fname = prefix + '_misaligned.csv'

  chrom_lens, mis = process_inputs(summary_fname, misalignments_fname, args['--down-sample'])
  plot_suffix = '.pdf' if args['pdf'] else '.png'

  if args['circle']:
    plot_fname = prefix + '_circle_plot' + plot_suffix
    draw_circle_plot(chrom_lens, mis)
  elif args['matrix']:
    plot_fname = prefix + '_matrix_plot' + plot_suffix
    draw_matrix_plot(chrom_lens, mis)

  plt.savefig(plot_fname)


def process_inputs(summary_fname, misalignments_fname, down_sample=None):
  """Read in the .json and .csv files produced by perfectbam and prepare the data structures we need for plotting

  :param summary_fname:        the _summary.json file
  :param misalignments_fname:  the _misalignments.csv file
  :param down_sample:          the factor by which to skip reads during loading
  :return: chrom_lens, misalignments as needed by the plotting functions
  """
  summary = json.load(open(summary_fname, 'r'))
  chrom_lens = [sh['LN'] for sh in summary['sequence_header']]
  with open(misalignments_fname, 'r') as csv_fp:
    csv_fp.readline()  # Get rid of the header
    if down_sample is None:
      misalignments = [[int(r) for r in [row[2], row[3], row[5], row[6]]] for row in csv.reader(csv_fp)]
    else:
      skip = int(down_sample)
      misalignments = [[int(r) for r in [row[2], row[3], row[5], row[6]]] for n, row in enumerate(csv.reader(csv_fp)) if n % skip == 0]

  return chrom_lens, misalignments


if __name__ == "__main__":
  cli()

