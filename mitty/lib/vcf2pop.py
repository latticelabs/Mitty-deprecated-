"""This module contains functions to parse a VCF file and return a master-list and chrom that can be used by the
rest of the system and to save the VCF as a genome db. It is used by `reads` to use VCF files directly as input"""
import gzip
import re
import io
import logging

import numpy as np

import mitty.lib.variants as vr


logger = logging.getLogger(__name__)


def vcf_to_pop(vcf_fname, pop_fname, sample_name=None, master_is_sample=False,
               genome_metadata=None,
               progress_callback=None, callback_interval=None):
  """This parses a VCF file into a Population object. Depending on how large your VCF is, you'll need a lot of memory.

  :param vcf_fname:
  :param pop_fname:
  :param sample_name: If None, the first sample will be taken. If not found, a runtime error is raised
  :param master_is_sample: If True we ignore any entires not part of the sample, resulting in a master list
                           that is identical and truncated to the sample
  :param genome_metadata: Pass this to override metadata in the VCF (or if it does not exist in the VCF)
  :param progress_callback: A function that will be called with the number of bytes read since the last call
  :param callback_interval: Lines to read before triggering callback
  :return: a Population object with a master list and one sample

  master list = all entries, unless sample_is_master = True
  sample = chosen sample

  Special cases:
  1. No samples (no GT column). Sample is same as master
  2. No sample name given. First sample in list is taken. If no samples, do what 1. requires
  """
  if progress_callback is not None:
    assert callback_interval is not None, "Callback interval must be set if progress callback function is set"

  with io.BufferedReader(gzip.open(vcf_fname, 'r')) if vcf_fname.endswith('gz') else open(vcf_fname, 'r') as fp:
    ftell = fp.raw.fileobj.tell if vcf_fname.endswith('gz') else fp.tell
    # Long story short, to get a reasonable indication of progress on a zipped file we need to know our actual
    # position in the compressed file, not uncompressed data stream. fp.tell() gets us the uncompressed position.
    # Since we handle both compressed and uncompressed files our ftell needs to be tailored to that
    _genome_metadata, gt_info_present, sample_column, sn = parse_header(fp, sample_name=sample_name)
    if genome_metadata is None: genome_metadata = _genome_metadata
    pop = vr.Population(fname=pop_fname, mode='w', genome_metadata=genome_metadata)
    for chrom, ml, svi in iter_vcf(
      fp,
      genome_metadata=genome_metadata,
      gt_info_present=gt_info_present,
      sample_column=sample_column,
      master_is_sample=master_is_sample,
      progress_callback=progress_callback,
      callback_interval=callback_interval,
      ftell=ftell
    ):
      pop.set_master_list(chrom, ml)
      pop.add_sample_chromosome(chrom, sample_name or 'anon', svi)
  return pop


def parse_header(fp, sample_name=None):
  """Given a file pointer, assume we are starting at the beginning of a VCF file and parse the whole header

  :param fp - file pointer
  :param sample_name
  :returns genome_metadata, gt_info_present, sample_column

  If the genotype field is absent -> gt_info_present=False
  If sample_name does not match, or no sample name is given, or genotype field is absent -> sample_column = None

  """
  def parse_genome_information(_line):
    ma = contig_re.findall(_line)[0].split(',')
    seq_id, seq_len, seq_md5 = 'None', 0, '0'
    for m in ma:
      val = m.split('=')
      if val[0] == 'ID': seq_id = val[1]
      elif val[0] == 'length': seq_len = int(val[1])
      elif val[0] == 'md5': seq_md5 = val[1]
    return {'seq_id': seq_id, 'seq_len': seq_len, 'seq_md5': seq_md5}

  def parse_column_header(_line, _sample_name):
    if not _line.startswith('#CHROM'):
      raise RuntimeError("Badly formatted VCF")

    cols = _line.split('\t')
    _gt_info_present = len(cols) > 8
    _sample_column, _sn = 8, 'anon'
    if _gt_info_present:
      if _sample_name is not None:
        s_c = [c for c, n in enumerate(cols) if n == _sample_name]
        if len(s_c) == 0:
          raise RuntimeError("No sample named {}".format(_sample_name))
        _sample_column, _sn = s_c[0], _sample_name
      else:
        _sample_column = 9
        _sn = cols[_sample_column]
    return _gt_info_present, _sample_column, _sn

  contig_re = re.compile(r"##contig=<(.*)>")
  genome_metadata, gt_info_present, sample_column, sn = [], False, -1, 'anon'
  for line in fp:
    if line[:8] == '##contig':
      genome_metadata.append(parse_genome_information(line.strip()))
    if line[:2] != '##':  # Done with the initial part of the header
      gt_info_present, sample_column, sn = parse_column_header(line.strip(), sample_name)
      break
  logger.debug('{} sequences found in VCF head ##contig information'.format(len(genome_metadata)))
  logger.debug('Sample name is "{}" ({})'.format(sn, ('col {}'.format(sample_column + 1) if sample_column > 8 else 'Inferred HOM')))
  return genome_metadata, gt_info_present, sample_column, sn


def iter_vcf(
  fp,
  genome_metadata=[],
  gt_info_present=False,
  sample_column=None,
  master_is_sample=False,
  progress_callback=None,
  callback_interval=None,
  ftell=None):
  """

  :param fp: pointer to VCF file
  :param genome_metadata: List of genome metadata as passed to Population
  :param gt_info_present: If this is False, we set ignore_genotype and sample_is_master (and ignore sample_column)
  :param sample_column: Which column of the VCF to get our sample GT from, if None, sample is same as master
  :param master_is_sample: If true, we only load the variants from the sample
  :param progress_callback: A function that will be called with the number of bytes read since the last call
  :param callback_interval: Lines to read before triggering callback
  :param ftell: function that returns actual position in file
  :return: An iterator over chrom, ml, svi.

  Special cases:
  1. No samples (no GT column). Sample is same as master
  2. No sample name given. First sample in list is taken. If no samples, do what 1. requires

  **NOTE** We assume the VCF is sorted by chrom and pos.
  """
  n2id = {v['seq_id']: i + 1 for i, v in enumerate(genome_metadata)}
  imprecise = ['<', '>', ':', '[', ']']
  sc = sample_column

  l_chrom, l_pos, l_stop, l_ref, l_alt, l_svi = -1, [], [], [], [], []
  if progress_callback: st_ln, f_pos = callback_interval, 0
  for line in fp:
    if progress_callback:
      st_ln -= 1
      if st_ln == 0:
        progress_callback(ftell() - f_pos)
        st_ln, f_pos = callback_interval, ftell()

    cols = line.split(None, sc + 1)
    this_chrom, pos, ref, _alts = n2id[cols[0]], int(cols[1]) - 1, cols[3], cols[4]

    if this_chrom != l_chrom:  # Time to flush!
      if l_pos:
        ml = vr.VariantList(l_pos, l_stop, l_ref, l_alt, np.ones(len(l_pos), dtype='f2'))
        ml.sorted = True
        yield l_chrom, ml, vr.l2ca(l_svi)

      l_chrom, l_pos, l_stop, l_ref, l_alt, l_svi = this_chrom, [], [], [], [], []

    if gt_info_present:  # Sample_column is guaranteed to have a valid value
      gt_col = cols[sc].split(':', 1)[0]
      if '.' in gt_col:  # No call
        h = [0, 0]
      elif '|' in gt_col:
        h = [int(_h) for _h in gt_col.split('|')]
      elif '/' in cols[sc]:
        h = [int(_h) for _h in gt_col.split('/')]
      else:
        h = [int(gt_col), int(gt_col)]  # These are calls on Y or male X Mitty treats them as homozygous for that chrom
    else:
      h = [1, 1]

    # Expand multi-allelic entries to bi-allelic ones and process
    for n, alt in enumerate(_alts.split(',')):
      # Skip any imprecise entries
      if any(imp in alt for imp in imprecise): continue
      if any(imp in ref for imp in imprecise): continue

      gt_match = False
      if n + 1 == h[1]:
        gt = 2 if n + 1 == h[0] else 1
        gt_match = True
      elif n + 1 == h[0]:
        gt = 0
        gt_match = True

      if not master_is_sample or gt_match:
        l_pos += [pos]
        l_stop += [pos + len(alt)]
        l_ref += [ref]
        l_alt += [alt]

      if gt_match:
        l_svi += [(len(l_pos) - 1, gt)]

  if progress_callback: progress_callback(ftell() - f_pos)
  if l_pos:
    ml = vr.VariantList(l_pos, l_stop, l_ref, l_alt, np.ones(len(l_pos), dtype='f2'))
    ml.sorted = True
    yield l_chrom, ml, vr.l2ca(l_svi)