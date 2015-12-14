"""This serves as the "meta-analysis" tool for aligner benchmarking, to collate the individual benchmarking
results and produce a summary html page."""
import csv
import json
from collections import OrderedDict
import logging

import click
import numpy as np

logger = logging.getLogger(__name__)


@click.command()
@click.argument('analysis_list', type=click.Path(exists=True))
@click.option('--out', type=click.Path(), help='Output file')
def cli(analysis_list, out):
  """Given a list of analyses as a CSV file load the files and metadata and create a summary page

  indel_anal, time, mis_plot, indel_plot, tool, graph, sample, read_set"""
  data = [row for row in csv.reader(open(analysis_list, 'r'), skipinitialspace = True)]
  header = data[0]
  assert header == ['indel_anal', 'time', 'mis_plot', 'indel_plot', 'tool', 'graph', 'sample', 'read_set']
  open(out, 'w').write(
    create_summary_page([[('tool', tool),
                          ('graph', graph),
                          ('sample', sample),
                          ('read_set', read_set),
                          ('indel_anal', summary_stats(load_indel_file(fn))),
                          ('time', parse_time_file(t_file)),
                          ('mis_plot', mis_plot),
                          ('indel_plot', indel_plot)]
                         for fn, t_file, mis_plot, indel_plot, tool, graph, sample, read_set in data[1:]])
  )


def load_indel_file(fn):
  """Load the json from the indel analysis

  :param fn: file name for the indel analysis
  :return: a dictionary with the json data converted into numpy arrays
  """
  data = json.load(open(fn, 'r'))
  data['indel_count'] = np.rec.fromrecords(data['indel_count'], dtype=[('x', 'int32'), ('total', 'uint32')])
  data['reads_within_feature'] = np.rec.fromrecords(data['reads_within_feature'], dtype=[('x', 'int32'), ('correct', 'uint32'), ('total', 'uint32')])
  return {'file': fn, 'data': data}


def summary_stats(f):
  """Simply tally up percentage corrects for Ref, SNP and Indel
    There are the following columns: DEL <= -20, -20 < DEL < 0, SNP, 0 < INS < 20, 20 <= INS
  """
  data = f['data']
  indel_len = data['indel_count']['x']
  _pc_raw = data['reads_within_feature']['correct'] / data['reads_within_feature']['total'].astype(float) * 100
  pc_raw = np.ma.masked_invalid(_pc_raw)

  ref_pc = data['fully_outside_features'][0] / float(data['fully_outside_features'][1]) * 100

  idx = np.nonzero(indel_len == 0)[0]
  snp_pc = pc_raw[idx][0]

  stats = {
    'ref': ref_pc,
    'snp': snp_pc
  }

  stats = []
  for k, r0, r1 in [('DEL <= -20', -1000, -20), ('-20 < DEL < 0', -20, -1),
                    ('0 < INS < 20', 0, 19), ('20 <= INS', 19, 1000)]:
    idx = np.nonzero((r0 < indel_len) & (indel_len <= r1))
    stats += [(k, pc_raw[idx].mean())]

  return OrderedDict(stats[0:2] + [('ref', ref_pc), ('snp', snp_pc)] + stats[2:])


def parse_time_file(t_file):
  """Expects the time file to have data in the format that is default for the unix time command

  time -p (ls -al) 2> time.txt

  real         0.00
  user         0.00
  sys          0.00

  We grab the 'real' entry
  """
  return float(open(t_file, 'r').readlines()[-3].strip().split()[1])


def parse_meta_data_file(m_file):
  """This is a json file that looks a bit like this:

  u'metadata': {u'INP:iset1': u'i1',
  u'INP:iset2': u'i3',
  u'bench_name': u'B1',
  u'bench_run': u'R5',
  u'file_extension': u'TXT',
  u'tool': u'tool1'}


  """
  return json.load(open(m_file, 'r')).get('metadata', {})


def create_summary_page(bench_stats):
  """Take a list of indel analysis and time stats and dump them to an html table.

  To start with, this is just a very flat html table:
  Each row is a tool run (unique combination of input files and tool).
  Eventually, the rows will be grouped by graph, sample, read_type and tool, right now, we just flatten the list
  There are N+3 columns.
  There are N columns for the indel categories described in the summary_stats function and the
  time, link to indel plot and link to mis-alignment plot

  :param bench_stats: Dictionary of values
  :return:
  """
  html = '<table border=1>'
  html += '<tr>'
  html += ''.join(['<th>{}</th>'.format(kk) for k, v in bench_stats[0]
                   for kk, _ in (v.items() if type(v) is OrderedDict else [(k, v)])])
  html += '</tr>\n<tr>'
  html += '</tr>\n<tr>'.join([''.join(['<td>{}</td>'.format(vv) for k, v in bs
                      for kk, vv in (v.items() if type(v) is OrderedDict else [(k, v)])]) for bs in bench_stats])
  html += '</tr>\n</table>'
  return html


if __name__ == '__main__':
  cli()