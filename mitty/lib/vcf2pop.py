"""This module contains functions to parse a VCF file and return a master-list and chrom that can be used by the
rest of the system and to save the VCF as a genome db. It is used by `reads` to use VCF files directly as input


Usage:
  vcf2pop --vcf=VCF --vdb=VDB --sample_name=SN

Options:
  --vcf=VCF    VCF file name
  --vdb=VDB    Genome database
  --sample_name=SN Name of sample
"""
import gzip
import re
import io

#import docopt
import numpy as np

import mitty.lib.variants as vr


def vcf_to_pop(vcf_fname, pop_fname, sample_name=None, master_is_sample=False):
  """This parses a VCF file into a Population object. Depending on how large your VCF is, you'll need a lot of memory.

  :param vcf_fname:
  :param pop_fname:
  :param sample_name: If None, the first sample will be taken. If not found, a runtime error is raised
  :param master_is_sample: If True we ignore any entires not part of the sample, resulting in a master list
                           that is identical and truncated to the sample
  :return: a Population object with a master list and one sample

  master list = all entries, unless sample_is_master = True
  sample = chosen sample

  Special cases:
  1. No samples (no GT column). Sample is same as master
  2. No sample name given. First sample in list is taken. If no samples, do what 1. requires
  """
  with io.BufferedReader(gzip.open(vcf_fname, 'r')) if vcf_fname.endswith('gz') else open(vcf_fname, 'r') as fp:
    genome_metadata, gt_info_present, sample_column = parse_header(fp, sample_name=sample_name)
    pop = vr.Population(fname=pop_fname, mode='w', genome_metadata=genome_metadata)
    for chrom, ml, svi in iter_vcf(
        fp,
        genome_metadata=genome_metadata,
        gt_info_present=gt_info_present,
        sample_column=sample_column,
        master_is_sample=master_is_sample):
      pop.set_master_list(chrom, ml)
      pop.add_sample_chromosome(chrom, sample_name, svi)
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
    if _gt_info_present:
      if _sample_name is not None:
        s_c = [c for c, n in enumerate(cols) if n == _sample_name]
        if len(s_c) == 0:
          raise RuntimeError("No sample named {}".format(_sample_name))
        _sample_column = s_c[0]
      else:
        _sample_column = 9
    else:
      _sample_column = None
    return _gt_info_present, _sample_column

  contig_re = re.compile(r"##contig=<(.*)>")
  genome_metadata, gt_info_present, sample_column = [], False, -1
  for line in fp:
    if line[:8] == '##contig':
      genome_metadata.append(parse_genome_information(line.strip()))
    if line[:2] != '##':  # Done with the initial part of the header
      gt_info_present, sample_column = parse_column_header(line.strip(), sample_name)
      break

  return genome_metadata, gt_info_present, sample_column


def iter_vcf(
  fp,
  genome_metadata=[],
  gt_info_present=False,
  sample_column=None,
  master_is_sample=False):
  """

  :param fp: pointer to VCF file
  :param genome_metadata: List of genome metadata as passed to Population
  :param gt_info_present: If this is False, we set ignore_genotype and sample_is_master (and ignore sample_column)
  :param sample_column: Which column of the VCF to get our sample GT from, if None, sample is same as master
  :param master_is_sample: If true, we only load the variants from the sample
  :return: An iterator over chrom, ml, svi.

  Special cases:
  1. No samples (no GT column). Sample is same as master
  2. No sample name given. First sample in list is taken. If no samples, do what 1. requires

  **NOTE** We assume the VCF is sorted by chrom and pos.
  """
  n2id = {v['seq_id']: i + 1 for i, v in enumerate(genome_metadata)}
  imprecise = ['<', '>', ':', '[', ']']

  l_chrom, l_pos, l_stop, l_ref, l_alt, l_svi = -1, [], [], [], [], []
  for line in fp:
    cols = line.split(None, sample_column)
    this_chrom, pos, ref, _alts = n2id[cols[0]], int(cols[1]) - 1, cols[3], cols[4]

    if this_chrom != l_chrom:  # Time to flush!
      if l_pos:
        ml = vr.VariantList(l_pos, l_stop, l_ref, l_alt, np.ones(len(l_pos), dtype='f2'))
        ml.sorted = True
        yield l_chrom, ml, vr.l2ca(l_svi)

      l_chrom, l_pos, l_stop, l_ref, l_alt, l_svi = this_chrom, [], [], [], [], []

    if gt_info_present:  # Sample_column is guaranteed to have a valid value
      h = [int(_h) for _h in (cols[-1].split('|') if '|' in cols[-1] else cols[-1].split('/'))]
    else:
      h = [0, 0]

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
        l_stop += [pos + len(alt)]  # Check this
        l_ref += [ref]
        l_alt += [alt]

      if gt_match:
        l_svi += [(len(l_pos) - 1, gt)]

  if l_pos:
    ml = vr.VariantList(l_pos, l_stop, l_ref, l_alt, np.ones(len(l_pos), dtype='f2'))
    ml.sorted = True
    yield l_chrom, ml, vr.l2ca(l_svi)



# def parse_vcf2(fname):
#   """This assumes a single sample VCF with only a GT field"""
#   opener = gzip.open if fname.endswith('gz') else open
#   with opener(fname) as fp:
#     header, seq_metadata = parse_header(fp)
#     master_lists = [vr.VariantList() for _ in header.keys()]
#     data = [[[], [], [], [], []] for _ in header.keys()]  # pos_a, stop_a, ref_a, alt_a, gt
#     zygosity = {'0|1': 1, '1|0': 0, '1|1': 2, '0/1': 1, '1/0': 0, '1/1': 2}
#     for line in fp:
#       cells = line.split()  # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s3
#       if cells[9] not in zygosity: continue
#       this_d = data[header[cells[0]]]
#       this_d[0].append(int(cells[1]) - 1)  #  VCF files are 1 indexed, we are, internally, 0 indexed
#       this_d[1].append(int(cells[1]) + len(cells[3]) - 1)
#       this_d[2].append(cells[3])
#       this_d[3].append(cells[4])
#       this_d[4].append(zygosity[cells[9]])
#
#   chroms = []
#   for n in range(len(header)):
#     this_d = data[n]
#     master_lists[n].add(this_d[0], this_d[1], this_d[2], this_d[3], [1.0] * len(this_d[0]))
#     master_lists[n].sort()
#     chroms.append(zip(range(len(this_d[4])), this_d[4]))
#   return master_lists, chroms, seq_metadata
#
#
# def vcf_to_pop2(vcf_fname, pop_fname, sample_name='s1', in_memory=False):
#   """Read a VCF file and store it as a Population structure
#   This assumes a single sample VCF with only a GT field"""
#   mls, chroms, genome_metadata = parse_vcf(vcf_fname)
#   pop = vr.Population(fname=pop_fname, mode='w', genome_metadata=genome_metadata, in_memory=in_memory)
#   for n, ml in enumerate(mls):
#     pop.set_master_list(n + 1, ml)
#     pop.add_sample_chromosome(n + 1, sample_name, np.array(chroms[n], dtype=[('index', 'i4'), ('gt', 'i1')]))
#   return pop


def cli():
  # args = docopt.docopt(__doc__)
  # vcf_to_pop(vcf_fname=args['--vcf'], pop_fname=args['--vdb'], sample_name=args['--sample_name'])
  pass


if __name__ == '__main__':
  cli()