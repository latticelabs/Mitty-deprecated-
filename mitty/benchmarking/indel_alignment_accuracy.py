"""Given a categorized reads file and

Usage:
  alindel <catreads> <pkl> (--vdb=VDB --sample_name=SN|--vcf=VCF) [--indel=INDEL]

Options:
  <catreads>   Categorized reads file
  <pkl>        name of output pkl file
  --vdb=VDB    Genome database
  --sample_name=SN Name of sample in genome database
  --vcf=VCF    VCF file name instead of genome database + sample_name
  --indel=INDEL  indel range [default: 100]"""
import cPickle

import docopt
import numpy as np

import mitty.lib.variants as vr
import mitty.lib.vcf2pop as vcf2pop
import mitty.benchmarking.creed as creed


def cli():
  args = docopt.docopt(__doc__)
  cat_fname=args['<catreads>']
  out_fname=args['<pkl>']
  max_indel = int(args['--indel'])
  if args['--vdb']:
    pop = vr.Population(fname=args['--vdb'])
    sample_name = args['--sample_name']
  else:
    sample_name = 's1'
    pop_fname = args['--vcf'] + '.pop.h5'
    pop = vcf2pop.vcf_to_pop(vcf_fname=args['--vcf'], pop_fname=pop_fname, sample_name=sample_name)
  cat_reads = creed.CategorizedReads(fname=cat_fname)
  ref_reads, cat_counts = process_sample(pop, cat_reads, sample_name, max_indel)
  cPickle.dump({'rr': ref_reads, 'cc': cat_counts}, open(out_fname, 'wb'))


def process_sample(pop, cat_reads, sample_name, max_indel):
  ref_reads = np.array([[0, 0], [0, 0]], dtype='i4')
  # rows -> one of pair is reference read, both of pair is reference read
  # cols -> correct, total
  cat_counts = None
  for ch in pop.get_chromosome_list():
    ml = pop.get_master_list(chrom=ch)
    vl = pop.get_sample_chromosome(chrom=ch, sample_name=sample_name)
    for cpy in [0, 1]:
      rl = cat_reads.get_data(chrom=ch, cpy=cpy)
      v_list = ml.variants[vl['index'][(vl['gt'] == cpy) | (vl['gt'] == 2)]]
      rc, vc0 = creed.bucket_list(rl['pos'], rl['stop'], rl['cat'], v_list['pos'], v_list['stop'])
      ref_reads += rc
      cat_counts = creed.categorize_read_counts_by_indel_length(v_list, vc0, cat_counts=cat_counts, max_indel=max_indel)
  return ref_reads, cat_counts


if __name__ == '__main__':
  cli()