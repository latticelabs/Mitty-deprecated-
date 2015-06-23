"""This module contains functions to parse a VCF file and return a master-list and chrom that can be used by the
rest of the system and to save the VCF as a genome db. It is used by `reads` to use VCF files directly as input"""
import gzip
import re

import mitty.lib.db as mdb
import mitty.lib.variants as vr


def parse_header(fp):
  """Given a file pointer, assume we are starting at the beginning of a VCF file and parse the whole header

  :param fp - file pointer
  :returns dict of chromosome names we can expect in col1 with values = to chromosome number
  """
  contig_re = re.compile(r"##contig=<(.*)>")
  header = {}
  seq_serial = 0
  seq_metadata = []
  for line in fp:
    if line[:2] != '##':
      break
    if line[:8] == '##contig':
      ma = contig_re.findall(line)[0].split(',')
      seq_id, seq_len, seq_md5 = 'None', 0, '0'
      for m in ma:
        val = m.split('=')
        if val[0] == 'ID': seq_id = val[1]
        elif val[0] == 'length': seq_len = int(val[1])
        elif val[0] == 'md5': seq_md5 = val[1]
      header[seq_id] = seq_serial
      seq_metadata.append({'seq_id': seq_id, 'seq_len': seq_len, 'seq_md5': seq_md5})
      seq_serial += 1
  return header, seq_metadata


def parse_vcf(fname):
  """This assumes a single sample VCF with only a GT field"""
  opener = gzip.open if fname.endswith('gz') else open
  with opener(fname) as fp:
    header, seq_metadata = parse_header(fp)
    master_lists = [vr.VariantList() for _ in header.keys()]
    data = [[[], [], [], [], []] for _ in header.keys()]  # pos_a, stop_a, ref_a, alt_a, gt
    zygosity = {'0|1': 1, '1|0': 0, '1|1': 2, '0/1': 1, '1/0': 0, '1/1': 2}
    for line in fp:
      cells = line.split()  # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s3
      if cells[9] not in zygosity: continue
      this_d = data[header[cells[0]]]
      this_d[0].append(int(cells[1]) - 1)  #  VCF files are 1 indexed, we are, internally, 0 indexed
      this_d[1].append(int(cells[1]) + len(cells[3]) - 1)
      this_d[2].append(cells[3])
      this_d[3].append(cells[4])
      this_d[4].append(zygosity[cells[9]])

  chroms = []
  for n in range(len(header)):
    this_d = data[n]
    master_lists[n].add(this_d[0], this_d[1], this_d[2], this_d[3], [1.0] * len(this_d[0]))
    master_lists[n].sort()
    chroms.append(zip(range(len(this_d[4])), this_d[4]))
  return master_lists, chroms, seq_metadata


def vcf_to_pop(vcf_fname, pop_fname):
  """Read a VCF file and store it as a Population structure
  This assumes a single sample VCF with only a GT field"""
  mls, chroms, seq_metadata = parse_vcf(vcf_fname)
  pop = vr.Population()


  conn = mdb.connect(db_name=db_fname)

  for n, meta in enumerate(seq_metadata):
    mdb.save_chromosome_metadata(conn, n + 1, **meta)
    if len(mls[n]):
      mdb.save_master_list(conn, n + 1, mls[n])
      if len(chroms[n]):
        mdb.save_sample(conn, 0, 0, n + 1, chroms[n])