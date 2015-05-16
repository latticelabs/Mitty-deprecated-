"""Plot mis-alignments data in multiple ways. Also serves as a library for querying a mis-alignments database

Commandline::

  Usage:
    plot_align circle <dbfile> [--fig=FIG] [--samples=S] [-i]
    plot_align matrix <dbfile> [--fig=FIG] [--samples=S]
    plot_align detailed <sfile> <dbfile>... [--fig=FIG] [--sample=S] [-i]

  Options:
    circle           Plot a circle plot
    matrix           Plot a 2D histogram
    detailed         Plot a detailed analysis panel
    <dbfile>         File containing mis alignment data as produced by perfectbam
    <sfile>          Genome score file
    --fig=FIG        Name of output figure file [default: align_fig.png]
    --samples=S      Samples to show [default: 1000]
"""
import os
import docopt
import json
import sqlite3 as sq

import matplotlib
orig_backend = matplotlib.get_backend()
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
import numpy as np

import mitty.util.perfectbam as pbam

import logging
logger = logging.getLogger(__name__)


# def process_inputs(summary_fname, misalignments_fname, down_sample=None):
#   """Read in the .json and .csv files produced by perfectbam and prepare the data structures we need for plotting
#
#   :param summary_fname:        the _summary.json file
#   :param misalignments_fname:  the _misalignments.csv file
#   :param down_sample:          the factor by which to skip reads during loading
#   :return: chrom_lens, mis_alignments as needed by the plotting functions
#   """
#   summary = json.load(open(summary_fname, 'r'))
#   chrom_lens = [sh['LN'] for sh in summary['sequence_header']]
#   with open(misalignments_fname, 'r') as csv_fp:
#     csv_fp.readline()  # Get rid of the header
#     if down_sample is None:
#       mis_alignments = [[int(r) for r in [row[2], row[3], row[5], row[6]]] for row in csv.reader(csv_fp)]
#     else:
#       skip = int(down_sample)
#       mis_alignments = [[int(r) for r in [row[2], row[3], row[5], row[6]]] for n, row in enumerate(csv.reader(csv_fp)) if n % skip == 0]
#
#   return chrom_lens, mis_alignments


class CircularPlot:
  """Encapsulates data and methods for the circular plot.
  In non-interactive mode we simply plot a sample of mis-aligned reads in the given axes.
  In interactive mode we use mouse-scroll to change the read-window size and we use hover to update which section of the
  genome we show

  When this plot is to be used as part of a larger, interactive plot, we should create this in non-interactive mode and
  call the functions from the governing controls
  """
  def __init__(self, conn, chrom=None, window=int(1e5), samples=1000, interactive=False, ax=None,
               radius=1, chrom_gap=0.05, chrom_thick=5, lw=1):
    """
    :param conn: database connection
    :param chrom: chromosome to start in. If None, we will plot all data
    :param window: window size of observable reads (ignored if chrom is None)
    :param samples: number of samples to show. If None, show all
    :param interactive: if True, create an interactive plot
    :param ax: axis object. Needs to be polar (e.g. ax = plt.subplot(111, polar=True))
    :param radius: radius of the plot
    :param chrom_gap: gap between chromosomes in radians
    :param chrom_thick: thickness of chromosome
    """
    self.conn = conn
    self.chrom = chrom
    self.pos = 0
    self.window = window
    self.samples = samples
    self.interactive = interactive
    self.radius = radius
    self.chrom_gap = chrom_gap
    self.chrom_thick = chrom_thick
    self.lw = lw

    self.summary = pbam.load_summary(self.conn)
    self.chrom_lens = [s['seq_len'] for s in self.summary]
    self.chrom_offsets = reduce(lambda x, y: x + [x[-1] + y], self.chrom_lens, [0])

    if ax is None:
      self.fig = plt.figure()
      self.ax = self.fig.add_subplot(111, polar=True)
    else:
      self.fig = ax.figure
      self.ax = ax

    self.current_patches = []  # We use this for the interactive plot
    self.plot_genome_as_a_circle()
    self.plot_reads()

    if interactive:
      self.cid = self.fig.canvas.mpl_connect('motion_notify_event', self.select_pos)
      plt.show()

  def select_pos(self, event):
    if event.inaxes != self.ax: return

    c_off, c_gap = self.chrom_offsets, self.chrom_gap
    total_len = sum(self.chrom_lens)
    radians_per_base = (2.0 * np.pi - len(self.chrom_lens) * c_gap) / total_len  # With allowance for chrom gaps

    theta = event.xdata
    for n, off in enumerate(c_off):
      if theta <= off * radians_per_base + n * c_gap:
        self.chrom = n
        self.pos = int((theta - (self.chrom - 1) * c_gap) / radians_per_base - c_off[self.chrom - 1])
        break
    self.update_plot()
    print self.chrom, self.pos

  def plot_genome_as_a_circle(self):
    """Plot the chromosomes on a circle."""
    total_len = sum(self.chrom_lens)
    radians_per_base = (2.0 * np.pi - len(self.chrom_lens) * self.chrom_gap) / total_len  # With allowance for chrom gaps

    xticks = []
    xticklabels = []
    delta_radian = 0.01
    start_radian = 0
    for ch_no, l in enumerate(self.chrom_lens):
      end_radian = start_radian + l * radians_per_base
      theta = np.arange(start_radian, end_radian, delta_radian)
      self.ax.plot(theta, [self.radius * 1.01] * theta.size, lw=self.chrom_thick)
      xticks.append((start_radian + end_radian)/2)
      xticklabels.append(str(ch_no + 1))
      start_radian = end_radian + self.chrom_gap

    plt.setp(self.ax.get_yticklabels(), visible=False)
    self.ax.grid(False)
    plt.setp(self.ax, xticks=xticks, xticklabels=xticklabels)
    self.ax.set_rmax(1.04 * self.radius)

  def plot_reads(self):
    c_off, c_gap, radius, lw, ax = self.chrom_offsets, self.chrom_gap, self.radius, self.lw, self.ax
    total_len = sum(self.chrom_lens)
    radians_per_base = (2.0 * np.pi - len(self.chrom_lens) * c_gap) / total_len  # With allowance for chrom gaps

    # http://matplotlib.org/users/path_tutorial.html
    codes = [
      Path.MOVETO,
      Path.CURVE4,
      Path.CURVE4,
      Path.CURVE4,
    ]

    patches_added = []
    for read in pbam.load_reads(self.conn, chrom=self.chrom,
                                start_pos=self.pos, stop_pos=self.pos + self.window, sample=self.samples):
      t0 = (c_off[read['correct_chrom']-1] + read['correct_pos']) * radians_per_base + (read['correct_chrom']-1) * c_gap
      t1 = (c_off[read['aligned_chrom']-1] + read['aligned_pos']) * radians_per_base + (read['aligned_chrom']-1) * c_gap
      this_radius = max(min(1.0, abs(t1 - t0) / np.pi), 0.05) * radius
      verts = [
        (t0, radius),  # P0
        (t0, radius - this_radius),  # P1
        (t1, radius - this_radius),  # P2
        (t1, radius),  # P3
      ]
      path = Path(verts, codes)
      patch = patches.PathPatch(path, facecolor='none', lw=lw)
      ax.add_patch(patch)
      patches_added.append(patch)
    if self.chrom is not None:
      t0 = (c_off[self.chrom-1] + self.pos) * radians_per_base + (self.chrom-1) * c_gap
      t1 = (c_off[self.chrom-1] + self.pos + self.window) * radians_per_base + (self.chrom-1) * c_gap
      verts = [
        (t0, 1.03 * radius),  # P0
        (t0, 1.03 * radius),  # P1
        (t1, 1.03 * radius),  # P2
        (t1, 1.03 * radius),  # P3
      ]
      path = Path(verts, codes)
      patch = patches.PathPatch(path, color='gray', lw=2)
      ax.add_patch(patch)
      patches_added.append(patch)

    self.current_patches = patches_added

  def update_plot(self):
    for p in self.current_patches:
      p.remove()
    self.plot_reads()
    plt.draw()


def cli():
  """Command line script entry point."""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
    return
  else:
    args = docopt.docopt(__doc__)

  plot_fname = args['--fig']

  if args['-i']:
    switch_to_interactive_backend()

  conn = pbam.connect_to_db(args['<dbfile>'][0])
  mappability_data = np.load(args['<sfile>']) if args['<sfile>'] else None
  if args['circle']:
    cp = CircularPlot(conn, samples=int(args['--samples']), interactive=args['-i'])

  plt.savefig(plot_fname)






def switch_to_interactive_backend():
  """By default we use Agg because we want to make static plots, but if we choose the interactive option we need to
  switch back to the original, interactive, backend."""
  if orig_backend not in matplotlib.rcsetup.interactive_bk:
    backend_to_use = matplotlib.rcsetup.interactive_bk[0]
    logger.warning('Original backend {:s} is not in the interactive backend list. Using first interactive backend found {:s}'.format(orig_backend, backend_to_use))
  else:
    backend_to_use = orig_backend
  plt.switch_backend(backend_to_use)


class InteractivePlot:
  def __init__(self, chrom_lens, mis):
    self.chrom = 1
    self.pos = 0
    self.window = 1e5
    self.chrom_lens = chrom_lens
    self.mis = mis
    self.fig = plt.figure()
    self.ax = self.fig.add_subplot(111, polar=True)
    self.current_patches = []
    plot_genome_as_a_circle(chrom_lens, radius=1, chrom_gap=0.05, chrom_thick=5, ax=self.ax)
    self.plot_data()

  def filter_reads(self):
    c, p, w = self.chrom, self.pos, self.window
    return [m for m in self.mis if m[0] == c and p < m[1] < p + w]

  def plot_data(self):
    for p in self.current_patches:
      p.remove()
    self.current_patches = plot_mis_alignments_on_a_circle(self.chrom_lens, self.filter_reads(), section=(self.chrom, self.pos, self.pos + self.window), radius=1, chrom_gap=0.05, lw=0.5, ax=self.ax)
    plt.draw()

  def inext(self, event):
    self.pos += self.window
    if self.pos >= self.chrom_lens[self.chrom - 1]:
      self.chrom += 1
      self.pos = 0
      if self.chrom > len(self.chrom_lens):
        self.chrom = 1
    self.plot_data()


def interactive_plot(args):
  prefix = args['<prefix>'][0]  # This is a surprising behavior from docopt
  summary_fname = prefix + '_summary.json'
  mis_alignments_fname = prefix + '_misaligned.csv'
  chrom_lens, mis = process_inputs(summary_fname, mis_alignments_fname, args['--down-sample'])
  switch_to_interactive_backend()
  iplot = InteractivePlot(chrom_lens, mis)
  cid = iplot.fig.canvas.mpl_connect('button_press_event', iplot.inext)
  plt.show()


def non_interactive_plots(args):
  prefix = args['<prefix>'][0]  # This is a surprising behavior from docopt
  summary_fname = prefix + '_summary.json'
  mis_alignments_fname = prefix + '_misaligned.csv'

  chrom_lens, mis = process_inputs(summary_fname, mis_alignments_fname, args['--down-sample'])
  plot_suffix = '.pdf' if args['pdf'] else '.png'

  if args['circle']:
    plot_fname = prefix + '_circle_plot' + plot_suffix
    draw_static_circle_plot(chrom_lens, mis)
  elif args['matrix']:
    plot_fname = prefix + '_matrix_plot' + plot_suffix
    draw_static_matrix_plot(chrom_lens, mis)
  elif args['detailed']:
    plot_fname = prefix + '_mappability_plot' + plot_suffix
    f = np.load(args['<sfile>'])
    draw_mappability_plot(chrom_lens, mis, f)
  else:
    logger.warning('Unhandled option')
    return

  plt.savefig(plot_fname)


def draw_static_circle_plot(chrom_lens, mis):
  """Draw a static circle plot with all the mis-alignments

  :param (list) chrom_lens: list of chromosome lengths
  :param (list) mis: list of tuples (correct_chrom, correct_pos, aligned_chrom, aligned_pos) Can be generator
  """
  ax = plt.subplot(111, polar=True)
  plot_genome_as_a_circle(chrom_lens, radius=1, chrom_gap=0.05, chrom_thick=5, ax=ax)
  plot_mis_alignments_on_a_circle(chrom_lens, mis, radius=1, chrom_gap=0.05, lw=0.5, ax=ax)
  #ax.set_rmax(1.04)


def plot_genome_as_a_circle(chrom_lens, radius=1, chrom_gap=0.001, chrom_thick=0.1, ax=None):
  """Plot the chromosomes on a circle. Will plot on the currently active axes.

  :param (list) chrom_lens: list of chromosome lengths
  :param (float) radius: radius of the plot
  :param (float) chrom_gap: gap between chromosomes in radians
  :param (float) chrom_thick: thickness of chromosome
  :param (object) ax: axis object. Needs to be polar (e.g. ax = plt.subplot(111, polar=True))
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
    ax.plot(theta, [radius * 1.01] * theta.size, lw=chrom_thick)
    xticks.append((start_radian + end_radian)/2)
    xticklabels.append(str(ch_no + 1))
    start_radian = end_radian + chrom_gap

  plt.setp(ax.get_yticklabels(), visible=False)
  ax.grid(False)
  plt.setp(ax, xticks=xticks, xticklabels=xticklabels)
  ax.set_rmax(1.04 * radius)


def plot_mis_alignments_on_a_circle(chrom_lens, misalignments, section=(), radius=1, chrom_gap=0.001, lw=2, ax=None):
  """Plot bezier curves indicating where the given misalignment reads originated and landed.

  :param chrom_lens: list of chromosome lengths
  :param misalignments: list of tuples (correct_chrom, correct_pos, aligned_chrom, aligned_pos) Can be generator
  :param section: tuple of the form (chrom, start, stop) indicating which section we are looking at. Draws a marker
                  omit to omit plotting of this marker
  :param radius: radius of the plot
  :param lw: thickness of drawn line
  :param ax: axis object. Needs to be polar (e.g. ax = plt.subplot(111, polar=True))
  :returns list of patches added. These can be `removed` as needed
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

  patches_added = []
  for m in misalignments:
    t0 = (sum(chrom_lens[:m[0]-1]) + m[1]) * radians_per_base + (m[0]-1) * chrom_gap
    t1 = (sum(chrom_lens[:m[2]-1]) + m[3]) * radians_per_base + (m[2]-1) * chrom_gap
    this_radius = max(min(1.0, abs(t1 - t0) / np.pi), 0.05) * radius
    verts = [
      (t0, radius),  # P0
      (t0, radius - this_radius),  # P1
      (t1, radius - this_radius),  # P2
      (t1, radius),  # P3
    ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', lw=lw)
    ax.add_patch(patch)
    patches_added.append(patch)
  if tuple:
    t0 = (sum(chrom_lens[:section[0]-1]) + section[1]) * radians_per_base + (section[0]-1) * chrom_gap
    t1 = (sum(chrom_lens[:section[0]-1]) + section[2]) * radians_per_base + (section[0]-1) * chrom_gap
    verts = [
      (t0, 1.03 * radius),  # P0
      (t0, 1.03 * radius),  # P1
      (t1, 1.03 * radius),  # P2
      (t1, 1.03 * radius),  # P3
    ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, color='gray', lw=2)
    ax.add_patch(patch)
    patches_added.append(patch)

    #plt.plot(t0, .99 * radius, '.')
  return patches_added


def draw_static_matrix_plot(chrom_lens, mis):
  """Draw a static matrix plot with all the mis-alignments shown

  :param (list) chrom_lens: list of chromosome lengths
  :param (list) mis: list of tuples (correct_chrom, correct_pos, aligned_chrom, aligned_pos) Can be generator
  """
  chrom_offsets = get_chrom_offsets(chrom_lens)
  ax = plt.gca()
  plt.setp(ax.get_xticklabels(), visible=False)
  plt.setp(ax.get_yticklabels(), visible=False)
  plt.setp(ax, xticks=chrom_offsets, yticks=chrom_offsets, ylabel='Correct position', xlabel='Aligned position', aspect=1)
  #plot_mis_alignment_matrix(chrom_offsets, mis, limits=[2.5e6, 3e6], histogram=False, ax=ax)
  plot_mis_alignment_matrix(chrom_offsets, mis, histogram=False, ax=ax)


def get_chrom_offsets(chrom_lens):
  return reduce(lambda x, y: x + [x[-1] + y], chrom_lens, [0])


def flatten_coordinates(chrom_offsets, mis):
  """Convert data in (chrom, pos) format to data in (pos') format, where pos' pretends all the chromosomes are laid out
  end to end

  :param (list) chrom_offsets: list of chromosome offsets
  :param (list) mis: list of tuples (correct_chrom, correct_pos, aligned_chrom, aligned_pos) Can be generator
  """
  return np.array([[chrom_offsets[r[0] - 1] + r[1], chrom_offsets[r[2] - 1] + r[3]] for r in mis])


def plot_mis_alignment_matrix(chrom_offsets, mis, limits=[], histogram=False, ax=None):
  """

  :param (list) chrom_lens: list of chromosome lengths
  :param (list) mis: list of tuples (correct_chrom, correct_pos, aligned_chrom, aligned_pos) Can be generator
  """
  flattened_coordinates = flatten_coordinates(chrom_offsets, mis)
  if not limits:
    limits = [0, chrom_offsets[-1]]

  if histogram:
    h, xe, ye = np.histogram2d(flattened_coordinates[:, 0], flattened_coordinates[:, 1], bins=100,
                               range=[limits, limits])
    #plt.imshow(h, origin="lower", extent=(0, chrom_offsets[-1], 0, chrom_offsets[-1]), cmap=plt.cm.gray_r, norm=LogNorm(), interpolation=None)
    hdl = plt.matshow(h, origin="lower", extent=(0, chrom_offsets[-1], 0, chrom_offsets[-1]), cmap=plt.cm.gray_r, norm=LogNorm(), interpolation=None)
  else:
    hdl = plt.scatter(flattened_coordinates[:, 1], flattened_coordinates[:, 0], s=1, marker='.')

  plt.setp(ax, xlim=[0, chrom_offsets[-1]], ylim=limits)  # Need to set this after setting ticks
  return hdl


def draw_mappability_plot(chrom_lens, misalignments, npz_file):
  """Plot locations of misaligned plots along-with k-mer scores computed and stored in a file

  """
  ax = plt.gca()
  draw_vertical_genome_mappability(chrom_lens, npz_file, ax)
  hh = bin_reads(chrom_lens, misalignments)
  draw_vertical_binned_reads(chrom_lens, hh, ax)


def bin_reads(chrom_lens, misalignments, source=True):
  """Make a histogram of read locations from the mis-alignment data

  :param chrom_lens: list of chromosome lengths
  :param misalignments: list of tuples (correct_chrom, correct_pos, aligned_chrom, aligned_pos)
  :param source: True if we want to bin according to the correct (source) position of the reads
                 False is we want to bin according to the aligned (destination) position of the reads
  :return:
  """
  chrom_pos = [[] for _ in range(len(chrom_lens))]
  bins = [{'range': (0, ch_len), 'bins': ch_len / 10000} for ch_len in chrom_lens]
  for r in misalignments:
    chrom_pos[r[0] - 1] += r[1:2]
  return [np.histogram(chrom_pos[n], density=True, **bins[n]) for n in range(len(chrom_lens))]


def draw_vertical_genome_mappability(chrom_lens, npz_file, ax):
  step = npz_file['arr_0'][0]
  chrom_max_len = max(chrom_lens)
  for n, ch_len in enumerate(chrom_lens):
    x = npz_file['arr_{:d}'.format(n + 1)] + 5000 * (n + 1)
    y = chrom_max_len - np.arange(x.shape[0], dtype=int) * step
    ax.plot(x, y, color=[0.7, 0.7, 0.7], lw=0.5)


def draw_vertical_binned_reads(chrom_lens, hh, ax):
  chrom_max_len = max(chrom_lens)
  for n, ch_len in enumerate(chrom_lens):
    x = -hh[n][0]/hh[n][0].max() * 2000 + 5000 * (n + 1)
    y = chrom_max_len - (hh[n][1][:-1] + hh[n][1][1:])/2.0
    ax.plot(x, y, color=[0.3, 0.3, 0.3], lw=0.25)


if __name__ == "__main__":
  cli()