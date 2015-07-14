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

import click
import docopt
import numpy as np

import mitty.lib.variants as vr
import mitty.lib.vcf2pop as vcf2pop
import mitty.benchmarking.creed as creed


@click.command()
@click.argument('cat_fname', type=click.Path(exists=True))
@click.argument('out_fname', type=click.Path())
@click.argument('vdb', type=click.Path(exists=True))  # , help='File name of genome database')
@click.argument('sample_name')  # , help='Name of sample')
@click.option('--sample_name2', help='Sample name indicating graph variants in database')
@click.option('--indel-range', help='Maximum base pair count of indels we analyze', type=int, default=50)
@click.option('--indel-range2', help='Maximum base pair count of indels we analyze', type=int, default=10)
def cli(cat_fname, out_fname, vdb, sample_name, sample_name2, indel_range, indel_range2):
  pop = vr.Population(vdb)
  cat_reads = creed.CategorizedReads(fname=cat_fname)
  ref_reads, cat_counts = process_sample(pop=pop, cat_reads=cat_reads, sample_name=sample_name, sample_name2=sample_name2, max_indel=indel_range, max_indel2=indel_range2)
  cPickle.dump({'rr': ref_reads, 'cc': cat_counts}, open(out_fname, 'wb'))

  # args = docopt.docopt(__doc__)
  # cat_fname=args['<catreads>']
  # out_fname=args['<pkl>']
  # max_indel = int(args['--indel'])
  # if args['--vdb']:
  #   pop = vr.Population(fname=args['--vdb'])
  #   sample_name = args['--sample_name']
  # else:
  #   sample_name = 's1'
  #   pop_fname = args['--vcf'] + '.pop.h5'
  #   pop = vcf2pop.vcf_to_pop(vcf_fname=args['--vcf'], pop_fname=pop_fname, sample_name=sample_name)
  # cat_reads = creed.CategorizedReads(fname=cat_fname)
  # ref_reads, cat_counts = process_sample(pop=pop, cat_reads=cat_reads, sample_name=sample_name, sammax_indel=max_indel)
  # cPickle.dump({'rr': ref_reads, 'cc': cat_counts}, open(out_fname, 'wb'))


def process_sample(pop, cat_reads, sample_name, sample_name2=None, max_indel=500, max_indel2=10, max_dist=30):
  ref_reads = np.array([[0, 0], [0, 0]], dtype='i4')
  # rows -> one of pair is reference read, both of pair is reference read
  # cols -> correct, total
  cat_counts = None
  cat_shared_counts = None
  for ch in pop.get_chromosome_list():
    ml = pop.get_master_list(chrom=ch)
    vl = pop.get_sample_chromosome(chrom=ch, sample_name=sample_name)
    vl2 = pop.get_sample_chromosome(chrom=ch, sample_name=sample_name2) if sample_name2 else None
    for cpy in [0, 1]:
      rl = cat_reads.get_data(chrom=ch, cpy=cpy)
      v_list = ml.variants[vl['index'][(vl['gt'] == cpy) | (vl['gt'] == 2)]]
      rc, vc0 = creed.bucket_list(rl['pos'], rl['stop'], rl['cat'], v_list['pos'], v_list['stop'])
      ref_reads += rc
      if sample_name2:
        v_list2 = ml.variants[vl2['index'][(vl2['gt'] == cpy) | (vl2['gt'] == 2)]]
        nearest_v2 = creed.find_nearest_variant(v_list, v_list2)
        cat_counts = creed.categorize_read_counts_by_indel_length_and_nearest_variant(
          v1=v_list, v_read_counts=vc0, nearest_v2=nearest_v2, cat_counts=cat_counts, max_v1_indel=max_indel,
          max_v2_indel=max_indel2, max_dist=max_dist)
      else:
        cat_counts = creed.categorize_read_counts_by_indel_length(variations=v_list, v_read_counts=vc0, cat_counts=cat_counts, max_indel=max_indel)
  return ref_reads, cat_counts


if __name__ == '__main__':
  cli()