"""Given a perfect BAM file from perfectbam use the information in the extended tags to work out alignment accuracy
parametrized by the sample indel lengths. Return information in a json file"""
from itertools import izip
import json

import click
import numpy as np
import pysam

import mitty.lib.variants as vr
import mitty.benchmarking.creed as creed


def prepare_features(pop, ch, sample_name=None):
  """Given a population database, a sample name and a chromosome, create a feature footprint vector to pass to
  creed.count_reads_under_features.

  :param pop: Population database
  :param ch: chromosome number (1, 2, 3 ...)
  :param sample_name: name of sample to take variant list from. Leave as None to take master list
  """
  sample_variant_list = [pop.get_variant_master_list(chrom=ch).variants] if sample_name is None else \
    pop.get_sample_variant_list_for_chromosome(chrom=ch, sample_name=sample_name)
  return [
    {'footprint': {'start': svl['pos'], 'stop': svl['stop']},
     'indel lengths': [len(a) - len(r) for a, r in izip(svl['alt'], svl['ref'])]}
    for svl in sample_variant_list]


def categorize_read_counts_by_indel_length(read_counts, indel_lengths, cat_read_counts=None, max_indel=100):
  assert max_indel > 0
  assert len(indel_lengths) == len(read_counts)
  if cat_read_counts is None:
    cat_read_counts = np.zeros(2 * max_indel + 1, dtype=[('x', 'int32'), ('correct', 'uint32'), ('total', 'uint32')])
    cat_read_counts['x'] = range(-max_indel, max_indel + 1)  # The range of indel lengths we are assuming.
  correct = cat_read_counts['correct']
  total = cat_read_counts['total']
  for rcc, rct, il in izip(read_counts['correct'], read_counts['total'], indel_lengths):
    if abs(il) > max_indel: continue
    correct[il + max_indel] += rcc
    total[il + max_indel] += rct
  return cat_read_counts


def categorize_indels_by_length(indel_lengths, cat_counts=None, max_indel=100):
  """Given an indel_length vector bin it by inde length and add to existing vector."""
  if cat_counts is None:
    cat_counts = np.zeros(2 * max_indel + 1, dtype=[('x', 'int32'), ('total', 'uint32')])
    cat_counts['x'] = range(-max_indel, max_indel + 1)  # The range of indel lengths we are assuming.

  indel_counts, _ = np.histogram(indel_lengths, bins=np.arange(-max_indel - 0.5, max_indel + 1.5),
                                 range=[-max_indel, max_indel])
  cat_counts['total'] += indel_counts.astype('uint32')
  return cat_counts


def categorize_data_from_one_chromosome(bam_fp, pop, ch, sample_name=None, cat_read_counts=None, max_indel=100):
  """For the given perfect BAM file categorize reads under the variants indicated and return the data binned by indel
  size

  :param bam_fp: input bam file as a pysam AlignmentFile class
  :param pop: Population database
  :param ch: chromosome number (1, 2, 3 ...)
  :param sample_name: name of sample to take variant list from. Leave as None to take master list
  :param max_indel: Longest indels to consider
  """
  if cat_read_counts is None: cat_read_counts = {'fully_outside_features': [0, 0],
                                                 'templates_within_feature_but_read_outside': None,
                                                 'reads_within_feature': None,
                                                 'indel_count': None}
  features = prepare_features(pop, ch, sample_name)
  f_chrom_id = bam_fp.header['SQ'][ch - 1]['SN']
  for cpy, f_v in enumerate(features):
    f_chrom_cpy = None if sample_name is None else cpy
    f_start = f_v['footprint']['start']
    f_stop = f_v['footprint']['stop']
    read_counts = creed.count_reads_under_features(bam_fp, f_chrom_id, f_start, f_stop, f_chrom_cpy=f_chrom_cpy)
    cat_read_counts['fully_outside_features'][0] += read_counts['fully_outside_features'][0]
    cat_read_counts['fully_outside_features'][1] += read_counts['fully_outside_features'][1]
    #for k in ['templates_within_feature_but_read_outside', 'reads_within_feature']:
    for k in ['reads_within_feature', 'templates_within_feature_but_read_outside']:
      cat_read_counts[k] = categorize_read_counts_by_indel_length(read_counts[k], f_v['indel lengths'],
                                                                  cat_read_counts=cat_read_counts[k],
                                                                  max_indel=max_indel)
    cat_read_counts['indel_count'] = categorize_indels_by_length(f_v['indel lengths'], cat_read_counts['indel_count'],
                                                                 max_indel=max_indel)

  return cat_read_counts


class NumpyJsonEncoder(json.JSONEncoder):
  def default(self, obj):
    if isinstance(obj, np.ndarray):
      return obj.tolist()
    return json.JSONEncoder.default(self, obj)


@click.command()
@click.argument('perbam', type=click.Path(exists=True))
@click.argument('vdb', type=click.Path(exists=True))  # , help='File name of genome database')
@click.argument('out-json', type=click.Path())
@click.option('--sample-name', help='Name of sample to compare against. Leave out to test against population')
@click.option('--indel-range', help='Maximum base pair count of indels we process', type=int, default=50)
def cli(perbam, out_json, vdb, sample_name, indel_range):
  """Compute alignment accuracy as a function of SNPs and indels in sample_name being covered by the reads in PERBAM.
  PERBAM should be from the perfectbam program. The result is stored as arrays in the OUT_JSON file"""
  bam_fp = pysam.AlignmentFile(perbam, 'rb')
  pop = vr.Population(vdb)
  cat_read_counts = None
  for ch in pop.get_chromosome_list():
    cat_read_counts = categorize_data_from_one_chromosome(
      bam_fp, pop, ch, sample_name=sample_name, cat_read_counts=cat_read_counts, max_indel=indel_range)

  json.dump(cat_read_counts, open(out_json, 'w'), cls=NumpyJsonEncoder)

if __name__ == '__main__':
  cli()