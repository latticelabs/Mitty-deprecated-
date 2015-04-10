"""Module that can run benchmarks and provides a base class for tool wrapping."""
import os
import json
from collections import OrderedDict
import logging

import pysam

import mitty.lib
from mitty.util import perfectbam

logger = logging.getLogger(__name__)


class DataSet(object):
  """Represents a data set."""
  def __init__(self, name, description):
    self.name = name
    self.description = description

def get_tool(name):
  tool = mitty.lib.load_benchmark_tool_wrapper(name)
  return tool.Tool


class Tool(object):
  """A wrapper round an aligner or data caller."""
  def __init__(self, name='No name', description='No description'):
    self.name = name
    self.description = description
    self.parameter_set = OrderedDict()

  def setup_files(self, file_set):
    """Run any steps needed to make the input files usable."""
    raise NotImplementedError

  def run(self, file_set, tool_p_set, out_prefix):
    """Run any pre-processing steps, execute the command line, run any post-processing steps.
    For aligner tools the return
    Return a dictionary
    with the output files."""
    raise NotImplementedError

  def add_parameter_set(self, **kwargs):
    raise NotImplementedError


# The information we need (about file sets, tool sets, parameters etc.) are easily saved in dictionaries, however the
# key names need to be specific. In order to make sure we use the right key names we use classes to wrap round the
# dictionary and provide a well defined interface


class BenchSet(object):
  """Everything needed for a benchmark"""
  def __init__(self):
    self.file_sets = OrderedDict()
    self.bench_p_sets = OrderedDict()
    self.tools = OrderedDict()
    self.exclude_list = []
    self.analyze_func = None

  def exclude_this(self, fid, tid, tpid):
    """Return True if the ids of this benchmark combination match in the exclude_list.
    :param (list) exclude_list: list of dictionaries with id combinations to be excluded.
                                If a value is None it matches to any
    :param fid: file set id
    :param tid: tool id
    :param tpid: tool parameter id
    :returns boolean
    """
    for excl in self.exclude_list:
      # if excl(.) is None it should match anything
      if (excl[0] or fid) == fid and (excl[1] or tid) == tid and (excl[2] or tpid) == tpid:
        return True
    return False

  def run_bench_marks(self, root_dir):
    """Run nested loops to execute all the benchmarks asked for. Create a multi dimensional dictionary that is saved to
    .json containing all the dimensions that we have in the benchmark space
    :param file_sets: dictionary of tuples  {id: (human name, file_set_dict) ... }
    :param bench_p_sets: dictionary of tuples  {id: (human name, bench_p_dict) ... }
    :param tools: dictionary of tuples {id: (human name, tool object) ... }
    :param tool_p_sets: dictionary of dictionary of tuples {tid: {id: (human name, tool_p_set) ...} ... }
                        tid needs to match with tool ids
    :param exclude_list: list of dictionaries with id combinations to be excluded.
    :param root_dir: where the output will be placed
    :param analyze_func: function that analyzes the results.
                         Typically:  analyze_bam or analyze_vcf
    """
    for fn, (fk, file_set) in enumerate(self.file_sets.iteritems()):
      logger.debug('Running file set {:s}'.format(fk))
      for tn, (tk, tool) in enumerate(self.tools.iteritems()):
        logger.debug('Running tool {:s}'.format(tk))
        tool.setup_files(file_set)
        for tpn, (tpk, _) in enumerate(tool.parameter_set.iteritems()):
          if self.exclude_this(fk, tk, tpk):
            continue
          out_prefix = os.path.join(root_dir, 'f{:d}_t{:d}_p{:d}'.format(fn, tn, tpn))
          if os.path.exists(out_prefix):
            raise RuntimeError('Output directory {:s} already exists'.format(out_prefix))
          os.makedirs(out_prefix)
          out_file_set = tool.run(file_set, tpk, out_prefix)
          for bpn, (bpk, bench_p_set) in enumerate(self.bench_p_sets.iteritems()):
            out_prefix2 = os.path.join(out_prefix, 'bps{:d}'.format(bpn))
            os.makedirs(out_prefix2)
            self.analyze_func(file_set=file_set, out_file_set=out_file_set, bench_p_set=bench_p_set, out_prefix=out_prefix2)

    bench_info = {
      'file_sets': self.file_sets.keys(),
      'bench_mark_sets': self.bench_p_sets.keys(),
      'tools': [{tn: self.tools[tn].parameter_set.keys()} for tn in self.tools.iterkeys()]
    }

    with open(os.path.join(root_dir, 'benchmark_info.json'), 'w') as fp:
      json.dump(bench_info, fp, indent=2)

    #generate_summary_jquery_page(root_dir, os.path.realpath('__file__'))


class FileSet(object):
  def __init__(self):
    pass


class ReferenceFileSet(FileSet):
  """Stores a reference file."""
  def __init__(self, fname):
    self.file = {'reference_file': os.path.join(os.path.abspath(__file__), fname)}



class AlignerFileSet(object):
  """Stores a dictionary with paths to the files."""
  def __init__(self, ):
    os.path.abspath()




class AlignerBenchSets(BenchSet):
  """Everything needed to benchmark an aligner."""
  def __init__(self):
    super(AlignerBenchSets, self).__init__()
    self.analyze_func = analyze_bam

  def add_file_set(self, name, reference_file, fastq, fastq_is_paired_interleaved):
    if name in self.file_sets:
      raise RuntimeError('This file set ({:s}) has already been used'.format(str(name)))
    self.file_sets[name] = {'reference_file': reference_file,
                            'fastq': fastq,
                            'fastq_is_paired_interleaved': fastq_is_paired_interleaved}

  def add_tool(self, tool):
    if tool.name in self.tools:
      raise RuntimeError('This tool name ({:s}) has already been used'.format(str(tool.name)))
    self.tools[tool.name] = tool

  def add_bench_parameter_set(self, name, window, block_size=1000000):
    if name in self.bench_p_sets:
      raise RuntimeError('This bench parameter name ({:s}) has already been used'.format(str(name)))
    self.bench_p_sets[name] = {'window': window, 'block_size': block_size}


def analyze_bam(out_file_set, bench_p_set, out_prefix, **kwargs):
  """
  :param file_set: ignored
  :param out_file_set: should have 'out_bam' filled in
  :param bench_p_set: should have 'window' and 'block_size'
  :param out_prefix: directory we should save our analysis files in"""

  assert 'out_bam' in out_file_set, '"out_file_set" dictionary returned by tool.run should have key "out_bam"'

  with pysam.Samfile(out_file_set['out_bam']) as bam_in_fp, \
       open(os.path.join(out_prefix, 'bench.csv'), 'w') as csv_out_fp, \
       open(os.path.join(out_prefix, 'bench.json'), 'w') as json_out_fp:
    perfectbam.main(bam_fp=bam_in_fp, csv_fp=csv_out_fp, json_fp=json_out_fp,
                  window=bench_p_set['window'], block_size=bench_p_set['block_size'])


def generate_summary_jquery_page(root_dir, script_dir):
  def summary_table(idx):
    html = '<table>'
    html += '<thead><tr><th rowspan="2">File sets</th>'
    for d in bench_info['tools']:
      tool_name = d.keys()[0]
      tool_param_names = d[tool_name]
      html += '<th colspan="{:d}">{:s}</th>'.format(len(tool_param_names), tool_name)
    html += '</tr>'
    html += '<tr>'
    for d in bench_info['tools']:
      for pn in d[d.keys()[0]]:
        html += '<th>{:s}</th>'.format(pn)
    html += '</tr></thead>'



    html += '</table>'
    return html

  includes = {
    'jquery-ui_min_css': open(os.path.join(script_dir, "jquery-ui-1.11.2/jquery-ui.min.css")).read(),
    'jquery_js': open(os.path.join(script_dir, "jquery-ui-1.11.2/external/jquery/jquery.js")).read(),
    'jqueryui_js': open(os.path.join(script_dir, "jquery-ui-1.11.2/jquery-ui.min.js")).read(),
    'jquery-ui_theme_min_css': open(os.path.join(script_dir, "jquery-ui-1.11.2/jquery-ui.theme.min.css")).read()
  }
  header_stuff = """
  <style media="screen" type="text/css">{jquery-ui_min_css}</style>
  <script>{jquery_js}</script>
  <script>{jqueryui_js}</script>
  <style media="screen" type="text/css">{jquery-ui_theme_min_css}</style>
  """.format(**includes)

  html_head = """<head>
  <meta charset="utf-8">
  <title>Benchmark results</title>""" + header_stuff + \
  """<script>
  $(function() {
    $( "#tabs" ).tabs();
  });
  </script>
  <style media="screen" type="text/css">
  table {
    color: #333;
    font-family: Helvetica, Arial, sans-serif;
    /*width: 640px;*/
    border-collapse: collapse;
    border-spacing: 0;
  }

  td, th {
    border: 1px solid; /* transparent; /* No more visible border */
    height: 30px;
    transition: all 0.3s;  /* Simple transition for hover effect */
  }

  th {
    background: #DFDFDF;  /* Darken header a bit */
    font-weight: bold;
  }

  td {
    background: #FAFAFA;
    text-align: center;
  }

  /* Cells in even rows (2,4,6...) are one color */
  tr:nth-child(even) td { background: #F1F1F1; }

  /* Cells in odd rows (1,3,5...) are another (excludes header cells)  */
  tr:nth-child(odd) td { background: #FEFEFE; }

  tr td:hover { background: #666; color: #FFF; } /* Hover cell effect! */

  /* Add border-radius to specific cells! */
  tr:first-child th:nth-child(2) {
    border-radius: 5px 0 0 0;
  }

  tr:first-child th:last-child {
    border-radius: 0 5px 0 0;
  }

  </style>
</head>"""

  with open(os.path.join(root_dir, 'benchmark_info.json'), 'r') as fp:
    bench_info = json.load(fp)
    html_body = '<body><div id="tabs"><ul>'
    for n, k in enumerate(bench_info['bench_mark_sets']):
      html_body += '<li><a href="#tabs-{:d}">{:s}</a></li>'.format(n, k)
    html_body += '</ul>'
    for n, k in enumerate(bench_info['bench_mark_sets']):
      html_body += '<div id="tabs-{:d}">'.format(n) + summary_table(n) + '</div>'
    html_body += '</div></body>'

  with open(os.path.join(root_dir, 'benchmark_info.html'), 'w') as fp:
    fp.write('<!doctype html><html lang="en">' + html_head + html_body + '</html>')
