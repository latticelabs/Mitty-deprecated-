"""Prepare a binned matrix of misalignments and plot it in different ways"""

import click
import pysam
import matplotlib
orig_backend = matplotlib.get_backend()
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
import numpy as np


def bam2mismat(bam_fp, bin_size=1000):
  """Create a matrix of binned mis-alignments"""
  def binnify(_pos, _bins):
    for n in range(1, len(_bins)):
      if _pos < _bins[n]:
        return n - 1
    return len(_bins) - 1 # Should not get here

  chrom_lens = [hdr['LN'] for hdr in bam_fp.header['SQ']]
  bins = [range(0, hdr['LN'], bin_size) + [hdr['LN']] for hdr in bam_fp.header['SQ']]
  # Rows = source (correct pos) Cols = destination (aligned pos)
  matrices = [[np.zeros(shape=(len(bins[j]) - 1, len(bins[i]) - 1), dtype='uint32') for i in range(len(bins))] for j in range(len(bins))]

  # TAG TYPE VALUE
  # XR  i    Aligned chromosome
  # XP  i    Aligned pos
  for r in bam_fp:
    c_chrom, c_pos, a_chrom, a_pos = r.reference_id, r.pos, r.get_tag('XR'), r.get_tag('XP')
    c_pos_binned, a_pos_binned = binnify(c_pos, bins[c_chrom]), binnify(a_pos, bins[a_chrom])
    matrices[c_chrom][a_chrom][c_pos_binned, a_pos_binned] += 1
  return chrom_lens, bins, matrices


def plot_genome_as_a_circle(ax, chrom_lens, chrom_gap=np.pi / 50, chrom_radius=1.0, chrom_thick=5, r_max=1.05):
  """Plot the chromosomes on a circle."""
  total_len = sum(chrom_lens)
  radians_per_base = (2.0 * np.pi - len(chrom_lens) * chrom_gap) / total_len  # With allowance for chrom gaps

  theta_stops, x_ticks, x_tick_labels = [], [], []
  delta_radian = 0.01
  start_radian = 0
  for ch_no, l in enumerate(chrom_lens):
    end_radian = start_radian + l * radians_per_base
    theta = np.arange(start_radian, end_radian, delta_radian)
    theta_stops.append((start_radian, end_radian))
    ax.plot(theta, [chrom_radius * 1.01] * theta.size, lw=chrom_thick)  # , color=[.3, .3, .3])
    x_ticks.append((start_radian + end_radian)/2)
    x_tick_labels.append(str(ch_no + 1))
    start_radian = end_radian + chrom_gap

  plt.setp(ax.get_yticklabels(), visible=False)
  ax.grid(False)
  plt.setp(ax, xticks=x_ticks, xticklabels=x_tick_labels)
  ax.set_rmax(r_max)
  return theta_stops


def plot_read_mis_alignments_on_a_circle(ax, chrom_lens, bins, matrices, theta_stops, chrom_radius=1.0, lw=0.01):
  bin_centers = [[(b0 + b1) / 2.0 for b0, b1 in zip(bb[:-1], bb[1:])] for bb in bins]

  # http://matplotlib.org/users/path_tutorial.html
  codes = [
    Path.MOVETO,
    Path.CURVE4,
    Path.CURVE4,
    Path.CURVE4,
  ]
  for i in range(len(bins)):
    for j in range(len(bins)):
      mat = matrices[i][j]
      for ii in range(mat.shape[0]):
        for jj in range(mat.shape[1]):
          if mat[ii, jj] == 0: continue

          t0 = theta_stops[i][0] + (bin_centers[i][ii] / float(chrom_lens[i])) * (theta_stops[i][1] - theta_stops[i][0])
          t1 = theta_stops[j][0] + (bin_centers[j][jj] / float(chrom_lens[j])) * (theta_stops[j][1] - theta_stops[j][0])

          this_radius = max(min(1.0, abs(t1 - t0) / np.pi), 0.05) * chrom_radius
          vertices = [
            (t0, chrom_radius),  # P0
            (t0, chrom_radius - this_radius),  # P1
            (t1, chrom_radius - this_radius),  # P2
            (t1, chrom_radius),  # P3
          ]
          path = Path(vertices, codes)
          patch = patches.PathPatch(path, facecolor='none', lw=lw * mat[ii, jj])
          ax.add_patch(patch)


def circle_plot(chrom_lens, bins, matrices):
  """Plot the confusion matrix as a circle plot."""
  fig = plt.figure()
  ax = fig.add_subplot(111, polar=True)
  theta_stops = plot_genome_as_a_circle(ax, chrom_lens)
  plot_read_mis_alignments_on_a_circle(ax, chrom_lens, bins, matrices, theta_stops, chrom_radius=1.0, lw=0.01)

@click.command()
@click.argument('badbam', type=click.Path(exists=True))
@click.option('--circle', type=click.Path(), help='Name of figure file for circle plot')
def cli(badbam, circle):
  """Prepare a binned matrix of mis-alignments and plot it in different ways"""
  chrom_lens, bins, matrices = bam2mismat(pysam.AlignmentFile(badbam, 'rb'), bin_size=10000)
  if circle is not None:
    circle_plot(chrom_lens, bins, matrices)
    plt.savefig(circle)

if __name__ == '__main__':
  cli()