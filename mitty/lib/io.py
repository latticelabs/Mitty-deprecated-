"""Some utilities related to loading/saving different file formats, compressing and indexing."""
import warnings
from os.path import splitext
import os
import glob
import gzip
from contextlib import contextmanager
import hashlib  # We decided to include md5 hashes of the sequences too
from itertools import izip

import pysam

import logging
logger = logging.getLogger(__name__)


MULTI_FASTA = 0
MULTI_DIR = 1


class Fasta:
  """This class handles loading of FASTA files.
  multi_fasta  -  a traditional gzipped fasta file storing multiple sequences, possibly with newlines
  multi_dir -  data split into multiple separate, unzipped, fasta files in one directory. Stored with no new-lines in
               the sequence. This allows us to speedily load individual sequences.
               This is useful if we have low memory (by setting persistent=False)
               or are prototyping and don't want to wait for ever to have the entire file to load.
               See the 'splitta' utility
  """
  def __init__(self, multi_fasta=None, multi_dir=None, persistent=True):
    """
    :param multi_fasta: fill out if input file is a single file
    :param multi_dir: fill out if input is in the form of multiple files in a directory numbered chr1.fa, chr2.fa etc.
    :param persistent: if True will keep sequences in memory after loading
    """
    assert not (multi_fasta is None and multi_dir is None), 'Need to specify either directory or file for reference. Check parameter file.'
    self.format = MULTI_DIR if multi_fasta is None else MULTI_FASTA
    self.multi_fasta = multi_fasta
    self.multi_dir = multi_dir

    self.sequences = {}  # This is a dict of dict of (seq, id, md5)
    self.persist = persistent

    if self.format == MULTI_DIR:
      self.seq_index = self.load_multi_dir_index()
      self._load_sequence_from_file = self.get_multi_dir
    else:
      self._load_sequence_from_file = self.get_multi_fasta
      _ = self[1]  # Just to load the sequences
      self.seq_index = [{'seq_id': self.sequences[n]['id'], 'seq_len': len(self.sequences[n]['seq']), 'seq_md5': self.sequences[n]['md5']}
                        for n in range(1, len(self.sequences) + 1)]
      if not self.persist:
        logger.warning('Persistence set to false for fa.gz file. Ignoring')

  def __getitem__(self, item):
    """This allows us to use Python's index notation to get sequences from the reference"""
    return self.sequences[item] if item in self.sequences else self._load_sequence_from_file(item)

  def load_multi_dir_index(self):
    """Load useful information about the genome from the index file.
    seqid, len and md5 sum
    """
    def dict_from_line(line):
      cells = line.split('\t')
      return {
        'seq_id': cells[0].strip(),
        'seq_len': cells[1].strip(),
        'seq_md5': cells[2].strip()
      }

    with open(glob.os.path.join(self.multi_dir, 'index.csv'), 'r') as fp:
      index = [dict_from_line(ln) for ln in fp]
    return index

  def get_multi_dir(self, item):
    """We get here because we don't have the sequence in memory"""
    fa_fname = glob.os.path.join(self.multi_dir, 'chr{:s}.fa'.format(str(item)))
    if not glob.os.path.exists(fa_fname):
      raise IOError('{:s} does not exist'.format(fa_fname))
    seq, sid = load_single_line_unzipped_fasta(fa_fname)
    ret_val = {'seq': seq, 'id': sid, 'md5': self.seq_index[item - 1]['seq_md5']}
    if self.persist:
      self.sequences[item] = ret_val
    return ret_val

  def get_multi_fasta(self, item):
    """We get here because we don't have the sequence in memory"""
    ref_seqs = load_generic_multi_fasta(self.multi_fasta)
    ret_val = {k: {'seq': v[0], 'id': v[1], 'md5': hashlib.md5(v[0]).hexdigest()} for k, v in ref_seqs.iteritems()}
    self.sequences = ret_val
    return ret_val[item]

  def __len__(self):
    return len(self.seq_index)

  def get(self, item, default=None):
    return self.__getitem__(item) or default

  def get_seq(self, chrom):
    return self[chrom]['seq']

  def get_seq_len(self, chrom):
    return len(self[chrom]['seq'])

  def get_seq_id(self, chrom):
    return self[chrom]['id']

  def get_seq_md5(self, chrom):
    return self[chrom]['md5']

  def get_seq_metadata(self):
    """Return a a list of (seq_id, seq_len, seq_md5) in same order as seen in fa.gz file"""
    return self.seq_index

  def __repr__(self):
    """Nice summary of what we have in the genome"""
    result = ''
    result += ('Multi fasta file: ' + self.multi_fasta) if self.multi_fasta else ('Multi dir: ' + self.multi_dir)
    result += '\n{:d} chromosomes\n'.format(len(self))
    for n in range(len(self)):
      result += '{:d} ({:s}) {:d} bases\n'.format(n + 1, self.get_seq_id(n + 1), self.get_seq_len(n + 1))
    return result

  def __iter__(self):
    for n in range(1, len(self) + 1):
      yield self[n]


def load_single_line_unzipped_fasta(fa_fname):
  """Expects a fasta file with only one sequence and only upper case letters - will read other files but the result
  is not sanitized in any way - newlines and repeat masks are left in.
  if as_numpy is set we will get the result as a numpy char array"""
  with open(fa_fname, 'r') as fasta_fp:
    seq_id = fasta_fp.readline()[1:-1]
    seq = fasta_fp.read()
  return seq, seq_id


def load_generic_multi_fasta(fa_fname):
  """Given a gzipped multi fa.gz file load it into a dictionary
  :param fa_fname: fasta.gz file with one or more fasta sequences
  :returns ref_seq: a dict of tuples of the form (seq, seq_id) with keys in the order they are found in the file

  Pure Python 2min 9s to load hg38
  Cythonized 2min 6s - since we are mostly in Python native functions, we are at speed limit
  """
  ref_seq = {}
  chr_no = 1
  with gzip.open(fa_fname, 'r') if fa_fname.endswith('gz') else open(fa_fname, 'r') as fp:
    l = fp.read()
  seq_strings = l.split('>')
  for seq_string in seq_strings:
    if seq_string == '':
      continue
    idx = seq_string.find('\n')
    if idx == -1:
      raise RuntimeError('Something wrong with the fasta file {:s}'.format(fa_fname))
    if idx == 0:
      continue  # Empty line, ignore
    seq_id = seq_string[:idx]
    seq = seq_string[idx:].replace('\n', '').upper()
    ref_seq[chr_no] = (seq, seq_id)
    chr_no += 1
  return ref_seq


# For now we concentrate on saving individual VCF files. Next version will have multi-vcf

@contextmanager
def vcf_for_writing(vcf_gz_name, sample_names, contig_info=[]):
  """Start a context with a handle to a vcf file. Write header before handing us the file pointer

  Given a master list for a given chromosome and a set of samples, save them to a VCF file.
  :param vcf_gz_name: Name of output vcf.gz file
  :param sample_names: list of sample names. If empty, write out master list with allele frequency
  :param contig_info: list of the form [(chrom, seq_id, seq_len, seq_md5) ... ]

  Example usage:
  with vcf_for_writing('test.vcf', ['a','b']) as fp:
    write_chromosomes_to_vcf(fp, chrom=1, chrom_list=[chrom_a, chrom_b], master_list=l)
  """
  if len(sample_names) > 1:
    raise NotImplementedError('Multiple sample VCF file saving is not supported in this version.')

  vcf_name, ext = splitext(vcf_gz_name)
  if ext != '.gz':
    vcf_name += '.vcf'
    vcf_gz_name = vcf_name + '.gz'

  header = "##fileformat=VCFv4.1\n"
  header += ''.join(['##contig=<ID={:s},length={:d},md5={:s}>\n'.format(ci[0].split(' ')[0], ci[1], ci[2]) for ci in contig_info])
  # When we write the contig id, we need to remove everything after the first space
  if len(sample_names):
    header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'  # We'll be writing samples
  else:
    header += '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'  # We'll be writing out master list
  header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  if len(sample_names): header += "\tFORMAT\t" + '\t'.join(sample_names)
  header += "\n"

  with open(vcf_name, 'w') as fp:
    fp.write(header)
    yield fp  # This where the writing happens

  compress_and_index_vcf(str(vcf_name), str(vcf_gz_name))


def compress_and_index_vcf(in_vcf_name, out_vcf_name):
  """Given an uncompressed, but sorted, vcf, compress and index it."""
  #bgzip -c sorted.vcf > sorted.vcf.gz
  #tabix sorted.vcf.gz
  logger.debug('Compressing and indexing {:s} to {:s}'.format(in_vcf_name, out_vcf_name))
  pysam.tabix_compress(in_vcf_name, out_vcf_name, force=True)
  pysam.tabix_index(out_vcf_name, force=True, preset='vcf')


def write_chromosomes_to_vcf(fp, seq_id='chr1', chrom_list=[], master_list=None):
  """Write out the chromosomes to as VCF lines

  :param fp: file pointer to context opened vcf file
  :param seq_id: sequence id, should match
  :param chrom_list: list of chromosome objects
  :param master_list: master list that the chromosome object indexes refer to
  """
  if len(chrom_list) > 1:
    raise NotImplementedError('Multiple sample VCF file saving is not supported in this version.')

  seq_id = seq_id.split(' ')[0]  # Only take contig_id upto up to the first space
  wr = fp.write
  gt_string = ['1|0', '0|1', '1|1']
  pos = master_list.variants['pos'] + 1  # VCF files are 1 indexed.
  ref = master_list.variants['ref']
  alt = master_list.variants['alt']
  maf = master_list.variants['p']

  if len(chrom_list) == 0:  # We want to write master list
    for p, r, a, f in izip(pos, ref, alt, maf):
      wr(seq_id + '\t' + str(p) + '\t.\t' + r + '\t' + a + '\t100\tPASS\tAF=' + str(f) + '\n')
  else:
    for idx, gt in chrom_list[0]:
      wr(seq_id + "\t" + str(pos[idx]) + "\t.\t" + ref[idx] + "\t" + alt[idx] + "\t100\tPASS\t.\tGT\t" + gt_string[gt] + "\n")


def load_vcf(rdr, sample_labels=None, max_rows=-1):
  """UNTESTED. SHOULD WORK

  Given a VCF reader from a particular chromosome, create master list and sample lists.
  :param rdr: VCF reader object (see below)
  :param sample_labels: list of sample names. If None, load all the samples in the file
  :returns l: master list of variants (dict)
  :returns samples: list of list of Genotypes

  fname = '/Users/kghose/Data/vcf1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz'
  vf = vcf.Reader(filename=fname)
  rdr = vf.fetch(chrom=22, start=0)
  """
  if sample_labels is None:
    sample_labels = rdr.samples

  l = {}
  samples = [[] for _ in sample_labels]
  for i, v in enumerate(rdr):
    if i == max_rows:  # Works for -1
      break
    for n, sn in enumerate(sample_labels):
      call = v.genotype(sn)
      if call.is_variant:
        g1 = int(call.gt_alleles[0])
        if g1 == 0:
          g1 = int(call.gt_alleles[1])
      else:
        continue
      g1 -= 1
      try:
        alt = v.ALT[g1].sequence
      except AttributeError:
        continue
      het = HOMOZYGOUS
      if call.gt_alleles[0] == '0':
        het = HET_01
      elif call.gt_alleles[1] == '0':
        het = HET_10
      samples[n] += [cgt(l, v.POS, stop=v.POS + len(v.REF), ref=v.REF, alt=alt, het=het)]
  return l, samples


def sort_and_index_bam(bamfile):
  """Do the filename gymnastics required to end up with a sorted, indexed, bam file."""
  # samtools sort adds a '.bam' to the end of the file name.
  os.rename(bamfile, 'temp.bam')
  pysam.sort('temp.bam', os.path.splitext(bamfile)[0])
  pysam.index(bamfile)
  os.remove('temp.bam')


def vcf2chrom(vcf_rdr):
  """Given a vcf reader corresponding to one chromosome, read in the data descriptions into our format. The result is
  sorted if the vcf file is sorted.
  """
  chrom = []  # deque()
  append = chrom.append
  for variant in vcf_rdr:
    alt = variant.ALT[0].sequence if variant.ALT[0] is not None else ''
    ref = variant.REF or ''
    start = variant.POS  # Note, we are in VCF coordinates!
    stop = variant.POS + len(ref)
    het = HOMOZYGOUS

    try:
      if variant.samples[0].gt_nums[0] == '0':
        het = HET_01
      if variant.samples[0].gt_nums[2] == '0':
        if het == HET_01:  # 0/0 means this does not exist in this sample
          het = ABSENT
          continue
        else:
          het = HET_10
    except IndexError:  # No genotype info, will assume homozygous
        pass

    append(new_variation(start, stop, ref, alt, het))

  return chrom


def parse_vcf(vcf_rdr, chrom_list):
  """Given a vcf reader load in all the chromosomes for all the samples."""
  g1 = {}
  for chrom in chrom_list:
    try:
      g1[chrom] = vcf2chrom(vcf_rdr.fetch(chrom, start=0))
    except (ValueError, KeyError):  # New version of pyvcf changed the error
      g1[chrom] = []  #deque()
  return g1


def vcf_save(g1, fp, sample_name='sample'):
  """Given a genome save it to a VCF file."""
  # Write header
  fp.write(
    "##fileformat=VCFv4.1\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{:s}\n".format(sample_name)
  )
  # Write lines
  wr = fp.write
  for chrom, variants in g1.iteritems():
    ch = str(chrom)
    for var in variants:
      # In the VCF file no REF or ALT is indicated by a .
      ref = var.vd.REF if var.vd.REF != '' else '.'
      alt = var.vd.ALT if var.vd.ALT != '' else '.'
      #  CHROM    POS   ID   REF   ALT   QUAL FILTER INFO FORMAT tsample
      wr(ch + "\t" + str(var.vd.POS) + "\t.\t" + ref + "\t" + alt + "\t100\tPASS\t.\tGT\t" + GT[var.het] + "\n")


def vcf_save_gz(g1, vcf_gz_name, sample_name='sample'):
  """Save .vcf, bgzip and index it. File name should have .gz at the end, but it's not a drama if doesnt. Sigh"""
  vcf_name, ext = splitext(vcf_gz_name)
  if ext != '.gz':  # Like I said, not a drama
    vcf_name += '.vcf'
    vcf_gz_name = vcf_name + '.gz'

  with open(vcf_name, 'w') as fp:
    vcf_save(g1, fp, sample_name=sample_name)

  compress_and_index_vcf(str(vcf_name), str(vcf_gz_name))
  # tabix can't understand unicode, needs bytes