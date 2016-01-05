"""This program chains together perfectbam, alindel, alindel_plot and misplot to do a complete analysis of an
aligned BAM file"""
import os
import subprocess

import click


@click.command()
@click.version_option()
@click.argument('inbam', type=click.Path(exists=True))
@click.argument('gdb', type=click.Path(exists=True))
@click.option('--cigar-errors', is_flag=True, help='CIGAR errors result in reads being classified as misaligned')
@click.option('--window', help='Size of tolerance window', default=0, type=int)
@click.option('--sample-name', help='Name of sample to compare against. Leave out to test against population')
@click.option('--indel-range', help='Maximum base pair count of indels we process', type=int, default=50)
@click.option('--indel-plot-smoothening-window', default=5, help='Size of median filter window to smooth plots')
@click.option('-p', is_flag=True, help='Show progress bar')
def cli(
  inbam, gdb,
  cigar_errors, window,
  sample_name, indel_range, indel_plot_smoothening_window, p):
  prog_bar = '-p' if p else ''

  prefix = os.path.splitext(os.path.basename(inbam))[0]

  bad_bam = prefix + '.bad.bam'
  per_bam = prefix + '.per.bam'

  # Perfectbam
  arguments = ['perfectbam', inbam, '-v', prog_bar, '--window', str(window)] \
              + (['--cigar-errors'] if cigar_errors else []) \
              + ['--bad-bam', bad_bam, '--per-bam', per_bam]
  subprocess.call(arguments)

  circle_plot = prefix + '.cir.pdf'
  matrix_plot = prefix + '.mat.pdf'

  # mis-plot
  arguments = ['misplot', '--circle', circle_plot, '--matrix', matrix_plot, bad_bam]
  subprocess.call(arguments)

  indel_json = prefix + '.indel.json'

  # alindel
  arguments = ['alindel', prog_bar, '--sample-name', sample_name, '--indel-range', str(indel_range), per_bam, gdb, indel_json]
  subprocess.call(arguments)

  alindel_plot = prefix + '.alindel.pdf'

  # alindel_plot
  arguments = ['alindel_plot', '-f', indel_json, '-o', alindel_plot, '-l', prefix,
               '--win', str(indel_plot_smoothening_window), '--indel-range', str(indel_range), '--title', 'indel: ' + prefix]
  subprocess.call(arguments)