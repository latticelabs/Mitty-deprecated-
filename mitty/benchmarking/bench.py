"""Module that can run benchmarks and provides a base class for tool wrapping."""
import os
from collections import OrderedDict

import pysam
from mitty.benchmarking import checkbam

import logging
logger = logging.getLogger(__name__)


class Tool(object):
  """A wrapper round an aligner or variant caller."""
  def __init__(self, name='No name'):
    self.name = name
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


class BenchResultMatrix:
  """A name indexed multi-dimensional table."""
  def __init__(self, bench_set):
    f_sets = bench_set.file_sets.iteritems()
    bp_sets = bench_set.bench_p_sets.iteritems()
    tools = bench_set.tools.iteritems()
    self.bench_mark_matrix = {
      bpk: {
        fsk: {
          tk: {
            tsk: tsv for tsk, tsv in tv.get_parameter_set().iteritems()
          } for tk, tv in tools
        } for fsk, _ in f_sets
      } for bpk, _ in bp_sets
    }

    # self.bench_mark_matrix = {bpk: [bench_p_sets[bpid][0],
    #                           {fid: [file_sets[fid][0],
    #                                  {tid: [tools[tid][0],
    #                                         {tpid: [tool_p_sets[tid][tpid][0], None]
    #                                          for tpk, tpv in tool_p_sets[tid])}]
    #                                  for tid in sorted(tools)}]
    #                           for fid in sorted(file_sets)}]
    #                    for bpid in sorted(bench_p_sets)}


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


def run_bench_marks2(file_sets, bench_p_sets, tools, tool_p_sets, exclude_list, root_dir, analyze_func):
  """Run nested loops to execute all the benchmarks asked for
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
  bench_mark_matrix = {bpid: [bench_p_sets[bpid][0],
                              {fid: [file_sets[fid][0],
                                     {tid: [tools[tid][0],
                                            {tpid: [tool_p_sets[tid][tpid][0], None]
                                             for tpid in sorted(tool_p_sets[tid])}]
                                     for tid in sorted(tools)}]
                              for fid in sorted(file_sets)}]
                       for bpid in sorted(bench_p_sets)}
  for fid, (fname, file_set) in file_sets.iteritems():
    logger.debug('Running file set {:s}'.format(fname))
    for tid, (tname, tool) in tools.iteritems():
      logger.debug('Running tool {:s}'.format(tname))
      tool.setup_files(file_set)
      for tpid, (tpname, tool_p_set) in tool_p_sets[tid].iteritems():
        if exclude_this(exclude_list, fid, tid, tpid):
          continue
        out_prefix = os.path.join(root_dir, 'f{:d}_t{:d}_p{:d}'.format(fid, tid, tpid))
        if os.path.exists(out_prefix):
          raise RuntimeError('Output directory {:s} already exists'.format(out_prefix))
        os.makedirs(out_prefix)
        out_file_set = tool.run(file_set, tool_p_set, out_prefix)
        for bpid, (bpname, bench_p_set) in bench_p_sets.iteritems():
          out_prefix2 = os.path.join(out_prefix, 'bps{:d}'.format(bpid))
          os.makedirs(out_prefix2)
          analyze_func(file_set=file_set, out_file_set=out_file_set, bench_p_set=bench_p_set, out_prefix=out_prefix2)
          bench_mark_matrix[bpid][1][fid][1][tid][1][tpid][1] = out_prefix2

  dimensions = {'file_sets': file_sets,
                'bench_p_sets': bench_p_sets,
                'tools': tools,
                'tool_p_sets': tool_p_sets}

  write_benchmark_metadata({'bench_mark_matrix': bench_mark_matrix, 'dimensions': dimensions}, root_dir)

html_wrapper = ("<!DOCTYPE html>\n"
                "<html lang=\"en\">\n"
                "  <head>\n"
                "    <title>Benchmarking results</title>\n"
                "  </head>\n"
                "  <body>\n"
                "  {:s}\n"
                "  </body>\n"
                "</html>\n")


def write_benchmark_metadata(bench_data, root_dir):
  """Save a .json file in root_dir that explains what benchmarks have been run. This information is helpful to create an
   index that links to all the benchmark results.



   """
  body = ""
  for bpid in sorted(bench_data['dimensions']['bench_p_sets']):
    body += '<table>'
    body += '<tr>'

    for tid in sorted(bench_data['dimensions']['tools']):
      body += '<th>'
      body += bench_data['dimensions']['tools'][tid][0]
      body += '</th>'
    body += '</tr>'

    for fid in sorted(bench_data['dimensions']['file_sets']):
      body += '<tr>'
      for tid in sorted(bench_data['dimensions']['tools']):
        for tpid in sorted(bench_data['dimensions']['tool_p_sets'][tid]):
          body += '<td>'
          body += bench_data['dimensions']['tool_p_sets'][tid][tpid][0]
          body += '</td>'
      body += '</tr>'
    body += '</table>'

  with open(os.path.join(root_dir, 'benchmark_summary.html'), 'w') as fp:
    fp.write(html_wrapper.format(body))


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
    checkbam.main(bam_fp=bam_in_fp, csv_fp=csv_out_fp, json_fp=json_out_fp,
                  window=bench_p_set['window'], block_size=bench_p_set['block_size'])  # TODO: Put these in bench mark parameters