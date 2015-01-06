"""Module that can run benchmarks and provides a base class for tool wrapping."""
import os
import pysam
from mitty.benchmarking import checkbam
import logging
logger = logging.getLogger(__name__)


class Tool:
  def __init__(self, _id=0, name='No name'):
    self.id = _id
    self.name = name

  def setup_files(self, file_set):
    """Run any steps needed to make the input files usable."""
    raise NotImplementedError

  def benchmark(self, file_set, tool_p_set, out_prefix):
    """Run any pre-processing steps, execute the command line, run any post-processing steps"""
    raise NotImplementedError


def exclude_this(exclude_list, fid, tid, tpid):
  """Return True if the ids of this benchmark combination match in the exclude_list.
  :param (list) exclude_list: list of dictionaries with id combinations to be excluded.
                              If a value is None it matches to any
  :param fid: file set id
  :param tid: tool id
  :param tpid: tool parameter id
  :returns boolean
  """
  for excl in exclude_list:
    # if excl(.) is None it should match anything
    if (excl[0] or fid) == fid and (excl[1] or tid) == tid and (excl[2] or tpid) == tpid:
      return True
  return False


def run_bench_marks(file_sets, tools, tool_p_sets, exclude_list, root_dir, analyze_func):
  """Run nested loops to execute all the benchmarks asked for
  :param file_sets: dictionary of tuples  {id: (human name, file_set_dict) ... }
  :param tools: dictionary of tuples {id: (human name, tool object) ... }
  :param tool_p_sets: dictionary of dictionary of tuples {tid: {id: (human name, tool_p_set) ...} ... }
                      tid needs to match with tool ids
  :param exclude_list: list of dictionaries with id combinations to be excluded.
  :param root_dir: where the output will be placed
  :param analyze_func: function that analyzes the results.
                       Typically:  analyze_bam or analyze_vcf
  """
  bench_mark_list = []
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
        out_file_set = tool.benchmark(file_set, tool_p_set, out_prefix)
        analyze_func(file_set, out_file_set, out_prefix)
        bench_mark_list += [{"id": (fid, tid, tpid),
                             "names": (fname, tname, tpname)}]
  create_result_page(bench_mark_list, root_dir)


def create_result_page(bench_mark_list, root_dir):
  pass


def analyze_bam(file_set, out_file_set, out_prefix):
  """
  :param file_set: ignored
  :param out_file_set: should have 'out_bam' filled in
  :param out_prefix: directory we should save our analysis files in"""

  with pysam.Samfile(out_file_set['out_bam']) as bam_in_fp, \
       open(os.path.join(out_prefix, 'bench.csv'), 'w') as csv_out_fp, \
       open(os.path.join(out_prefix, 'bench.json'), 'w') as json_out_fp:
    checkbam.main(bam_fp=bam_in_fp, csv_fp=csv_out_fp, json_fp=json_out_fp,
                  window=0, block_size=int(1e6))  # TODO: Put these in bench mark parameters