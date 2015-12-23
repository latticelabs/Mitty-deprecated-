import numpy as np

import mitty.lib.vcf2pop as vcf2pop
import mitty.lib.io as mio
import mitty.lib.variants as vr
import mitty.tests

import os
import tempfile

genome_metadata = [
  {'seq_id': 'NC_010127.1', 'seq_len': 422616, 'seq_md5': 'fe4be2f3bc5a7754085ceaa39a5c0414'},
  {'seq_id': 'NC_010128.1', 'seq_len': 457013, 'seq_md5': '99880025dcbcba5dbf72f437092903c3'},
  {'seq_id': 'NC_010129.1', 'seq_len': 481791, 'seq_md5': 'a3a5142f08b313f645cd5e972f5f3397'},
  {'seq_id': 'NC_010130.1', 'seq_len': 513455, 'seq_md5': 'ec8ff24820287d35c2b615fbb0df721c'},
]

variant_data = [
  {'pos': [1, 100, 200], 'stop': [2, 101, 201], 'ref': ['A', 'C', 'G'], 'alt': ['G', 'T', 'C'], 'p': [0.5, 0.5, 0.5]},
  {'pos': [1, 100, 200], 'stop': [2, 101, 201], 'ref': ['A', 'C', 'G'], 'alt': ['G', 'T', 'C'], 'p': [0.5, 0.5, 0.5]},
]

genotype_data = [
  [(0, 2), (2, 0)],
  [(1, 2), (2, 0)],
]


def round_trip_test():
  """vcf <-> mitty database round trip"""
  master_lists = {n + 1: vr.VariantList(vd['pos'], vd['stop'], vd['ref'], vd['alt'], vd['p']) for n, vd in enumerate(variant_data)}
  for k, v in master_lists.iteritems(): v.sort()

  pop = vr.Population(mode='w', genome_metadata=genome_metadata, in_memory=True)
  for k, v in master_lists.iteritems():
    pop.set_master_list(chrom=k, master_list=v)

  for n in [0, 1]:
    pop.add_sample_chromosome(n + 1, 'brown_fox', np.array(genotype_data[n], dtype=[('index', 'i4'), ('gt', 'i1')]))

  _, vcf_temp = tempfile.mkstemp(dir=mitty.tests.data_dir, suffix='.vcf.gz')
  _, h5_temp = tempfile.mkstemp(dir=mitty.tests.data_dir, suffix='.h5')

  mio.write_single_sample_to_vcf(pop, out_fname=vcf_temp, sample_name='brown_fox')

  pop2 = vcf2pop.vcf_to_pop(vcf_temp, h5_temp, sample_name='brown_fox')

  # for k, v in master_lists.iteritems():
  #   assert_array_equal(pop.get_sample_variant_index_for_chromosome(k, 'brown_fox'), pop2.get_sample_variant_index_for_chromosome(k, 'brown_fox'))

  for n in [0, 1]:  # Chromosomes
    for v1, v2 in zip(pop2.get_sample_variant_index_for_chromosome(n + 1, 'brown_fox'), genotype_data[n]):
      assert v1[1] == v2[1]  # Genotypes match
      for k in ['pos', 'stop', 'ref', 'alt']:
        assert pop2.get_variant_master_list(n + 1).variants[v1[0]][k] == master_lists[n + 1].variants[v2[0]][k]  # Variant data match

  os.remove(vcf_temp)