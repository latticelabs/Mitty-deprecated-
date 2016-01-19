import numpy as np
from numpy.testing import assert_array_equal

import mitty.lib.vcf2pop as vcf2pop
import mitty.lib.mio as mio
import mitty.lib.variants as vr
import mitty.tests

import os
import tempfile


def round_trip_test():
  """vcf <-> mitty database round trip"""
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


def vcf_reader_test1():
  """VCF with no GT data"""
  _vcf = """##fileformat=VCFv4.1
##contig=<ID=NC_010142.1,length=908485,md5=9a28f270df93bb4ac0764676de1866b3>
##contig=<ID=NC_010143.1,length=1232258,md5=ab882206d71bc36051f437e66246da6b>
##contig=<ID=NC_010144.1,length=1253087,md5=ab11fdfc260a2b78fdb845d89c7a89f2>
##contig=<ID=NC_010145.1,length=1282939,md5=b3c4b1a7b3671e2e8d4f4b1d2b599c44>
##contig=<ID=NC_010146.1,length=1621617,md5=3dbe62009f563fd1a6e3eadc15617e5c>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
NC_010142.1\t100\t.\tA\tT\t50\tPASS\tRS=672601345;RSPOS=1014319;VP=0x050060001205000002110200;GENEINFO=ISG15:9636;dbSNPBuildID=142;SAO=1;SSR=0;WGT=1;VC=DIV;PM;NSF;REF;ASP;LSD;OM;CLNHGVS=NC_000001.11:g.1014319dupG;CLNALLE=1;CLNSRC=OMIM_Allelic_Variant;CLNORIGIN=1;CLNSRCID=147571.0002;CLNSIG=5;CLNDSDB=MedGen:OMIM:Orphanet;CLNDSDBID=CN221808:616126:ORPHA319563;CLNDBN=Immunodeficiency_38;CLNREVSTAT=no_assertion_criteria_provided;CLNACC=RCV000148989.5
NC_010142.1\t120\t.\tA\tACT\t50\tPASS\tRS=672601345;RSPOS=1014319;VP=0x050060001205000002110200;GENEINFO=ISG15:9636;dbSNPBuildID=142;SAO=1;SSR=0;WGT=1;VC=DIV;PM;NSF;REF;ASP;LSD;OM;CLNHGVS=NC_000001.11:g.1014319dupG;CLNALLE=1;CLNSRC=OMIM_Allelic_Variant;CLNORIGIN=1;CLNSRCID=147571.0002;CLNSIG=5;CLNDSDB=MedGen:OMIM:Orphanet;CLNDSDBID=CN221808:616126:ORPHA319563;CLNDBN=Immunodeficiency_38;CLNREVSTAT=no_assertion_criteria_provided;CLNACC=RCV000148989.5
NC_010142.1\t140\t.\tACT\tA\t50\tPASS\tRS=672601345;RSPOS=1014319;VP=0x050060001205000002110200;GENEINFO=ISG15:9636;dbSNPBuildID=142;SAO=1;SSR=0;WGT=1;VC=DIV;PM;NSF;REF;ASP;LSD;OM;CLNHGVS=NC_000001.11:g.1014319dupG;CLNALLE=1;CLNSRC=OMIM_Allelic_Variant;CLNORIGIN=1;CLNSRCID=147571.0002;CLNSIG=5;CLNDSDB=MedGen:OMIM:Orphanet;CLNDSDBID=CN221808:616126:ORPHA319563;CLNDBN=Immunodeficiency_38;CLNREVSTAT=no_assertion_criteria_provided;CLNACC=RCV000148989.5
NC_010146.1\t100\t.\tA\tT\t50\tPASS\tRS=672601345;RSPOS=1014319;VP=0x050060001205000002110200;GENEINFO=ISG15:9636;dbSNPBuildID=142;SAO=1;SSR=0;WGT=1;VC=DIV;PM;NSF;REF;ASP;LSD;OM;CLNHGVS=NC_000001.11:g.1014319dupG;CLNALLE=1;CLNSRC=OMIM_Allelic_Variant;CLNORIGIN=1;CLNSRCID=147571.0002;CLNSIG=5;CLNDSDB=MedGen:OMIM:Orphanet;CLNDSDBID=CN221808:616126:ORPHA319563;CLNDBN=Immunodeficiency_38;CLNREVSTAT=no_assertion_criteria_provided;CLNACC=RCV000148989.5
"""
  vcf_name = os.path.join(mitty.tests.data_dir, 'vcf1.vcf')
  db_name = os.path.join(mitty.tests.data_dir, 't1.h5')
  open(vcf_name, 'w').write(_vcf)
  p = vcf2pop.vcf_to_pop(vcf_fname=vcf_name, pop_fname=db_name)

  assert p.get_chromosome_list() == [1, 2, 3, 4, 5], p.get_chromosome_list()
  assert p.get_sample_names() == ['anon']

  vl = p.get_sample_variant_list_for_chromosome(1,'anon')
  assert_array_equal(vl[0], vl[1], vl)
  assert vl[0][0]['pos'] == 99

  vl = p.get_sample_variant_list_for_chromosome(2, 'anon')
  assert vl[0].size == 0

  vl = p.get_sample_variant_list_for_chromosome(5,'anon')
  assert_array_equal(vl[0], vl[1], vl)


def vcf_reader_test2():
  """VCF with only GT data"""
  _vcf = """##fileformat=VCFv4.1
##contig=<ID=NC_010142.1,length=908485,md5=9a28f270df93bb4ac0764676de1866b3>
##contig=<ID=NC_010143.1,length=1232258,md5=ab882206d71bc36051f437e66246da6b>
##contig=<ID=NC_010144.1,length=1253087,md5=ab11fdfc260a2b78fdb845d89c7a89f2>
##contig=<ID=NC_010145.1,length=1282939,md5=b3c4b1a7b3671e2e8d4f4b1d2b599c44>
##contig=<ID=NC_010146.1,length=1621617,md5=3dbe62009f563fd1a6e3eadc15617e5c>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1
NC_010142.1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0|1
NC_010142.1\t120\t.\tA\tACT\t50\tPASS\t.\tGT\t1|0
NC_010142.1\t140\t.\tACT\tA\t50\tPASS\t.\tGT\t1|1
NC_010146.1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t.
"""
  vcf_name = os.path.join(mitty.tests.data_dir, 'vcf2.vcf')
  db_name = os.path.join(mitty.tests.data_dir, 't2.h5')
  open(vcf_name, 'w').write(_vcf)
  p = vcf2pop.vcf_to_pop(vcf_fname=vcf_name, pop_fname=db_name)

  vl = p.get_sample_variant_list_for_chromosome(1, 's1')
  assert vl[0][0]['pos'] == 119
  assert vl[0][0]['stop'] == 120
  assert vl[1][1]['pos'] == 139
  assert vl[1][1]['stop'] == 142

  assert vl[0].size == 2
  assert vl[1].size == 2


def vcf_reader_test3():
  """VCF with GT data and some other fields"""
  _vcf = """##fileformat=VCFv4.1
##contig=<ID=NC_010142.1,length=908485,md5=9a28f270df93bb4ac0764676de1866b3>
##contig=<ID=NC_010143.1,length=1232258,md5=ab882206d71bc36051f437e66246da6b>
##contig=<ID=NC_010144.1,length=1253087,md5=ab11fdfc260a2b78fdb845d89c7a89f2>
##contig=<ID=NC_010145.1,length=1282939,md5=b3c4b1a7b3671e2e8d4f4b1d2b599c44>
##contig=<ID=NC_010146.1,length=1621617,md5=3dbe62009f563fd1a6e3eadc15617e5c>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GL,Number=3,Type=Float,Description="Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1
NC_010142.1\t100\t.\tA\tT\t50\tPASS\t.\tGT:GL\t0|1:41,0,57
NC_010142.1\t120\t.\tA\tACT\t50\tPASS\t.\tGT:GL\t1|0:41,0,57
NC_010142.1\t140\t.\tACT\tA\t50\tPASS\t.\tGT:GL\t1|1:41,0,57
NC_010146.1\t100\t.\tA\tT\t50\tPASS\t.\tGT:GL\t.:41,0,57
"""
  vcf_name = os.path.join(mitty.tests.data_dir, 'vcf2.vcf')
  db_name = os.path.join(mitty.tests.data_dir, 't2.h5')
  open(vcf_name, 'w').write(_vcf)
  p = vcf2pop.vcf_to_pop(vcf_fname=vcf_name, pop_fname=db_name)

  vl = p.get_sample_variant_list_for_chromosome(1, 's1')
  assert vl[0][0]['pos'] == 119
  assert vl[0][0]['stop'] == 120
  assert vl[1][1]['pos'] == 139
  assert vl[1][1]['stop'] == 142

  assert vl[0].size == 2
  assert vl[1].size == 2


def vcf_reader_test4():
  """VCF with no genome metadata"""
  _vcf = """##fileformat=VCFv4.1
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GL,Number=3,Type=Float,Description="Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1
NC_010142.1\t100\t.\tA\tT\t50\tPASS\t.\tGT:GL\t0|1:41,0,57
NC_010142.1\t120\t.\tA\tACT\t50\tPASS\t.\tGT:GL\t1|0:41,0,57
NC_010142.1\t140\t.\tACT\tA\t50\tPASS\t.\tGT:GL\t1|1:41,0,57
NC_010146.1\t100\t.\tA\tT\t50\tPASS\t.\tGT:GL\t.:41,0,57
"""
  vcf_name = os.path.join(mitty.tests.data_dir, 'vcf2.vcf')
  db_name = os.path.join(mitty.tests.data_dir, 't2.h5')
  open(vcf_name, 'w').write(_vcf)

  genome_metadata = [
    {'seq_id': 'NC_010142.1', 'seq_len': 908485, 'seq_md5': '9a28f270df93bb4ac0764676de1866b3'},
    {'seq_id': 'NC_010143.1', 'seq_len': 1232258, 'seq_md5': 'ab882206d71bc36051f437e66246da6b'},
    {'seq_id': 'NC_010144.1', 'seq_len': 1253087, 'seq_md5': 'ab11fdfc260a2b78fdb845d89c7a89f2'},
    {'seq_id': 'NC_010145.1', 'seq_len': 1282939, 'seq_md5': 'b3c4b1a7b3671e2e8d4f4b1d2b599c44'},
    {'seq_id': 'NC_010146.1', 'seq_len': 1621617, 'seq_md5': '3dbe62009f563fd1a6e3eadc15617e5c'}
  ]

  p = vcf2pop.vcf_to_pop(vcf_fname=vcf_name, pop_fname=db_name, genome_metadata=genome_metadata)

  vl = p.get_sample_variant_list_for_chromosome(1, 's1')
  assert vl[0][0]['pos'] == 119
  assert vl[0][0]['stop'] == 120
  assert vl[1][1]['pos'] == 139
  assert vl[1][1]['stop'] == 142

  assert vl[0].size == 2
  assert vl[1].size == 2


def vcf_reader_test5():
  """VCF with multiple samples"""
  _vcf = """##fileformat=VCFv4.1
##contig=<ID=NC_010142.1,length=908485,md5=9a28f270df93bb4ac0764676de1866b3>
##contig=<ID=NC_010143.1,length=1232258,md5=ab882206d71bc36051f437e66246da6b>
##contig=<ID=NC_010144.1,length=1253087,md5=ab11fdfc260a2b78fdb845d89c7a89f2>
##contig=<ID=NC_010145.1,length=1282939,md5=b3c4b1a7b3671e2e8d4f4b1d2b599c44>
##contig=<ID=NC_010146.1,length=1621617,md5=3dbe62009f563fd1a6e3eadc15617e5c>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2
NC_010142.1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0|1\t1|0
NC_010142.1\t120\t.\tA\tACT\t50\tPASS\t.\tGT\t1|0\t0|1
NC_010142.1\t140\t.\tACT\tA\t50\tPASS\t.\tGT\t1|1\t0|0
NC_010146.1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t.\t1|1
"""
  vcf_name = os.path.join(mitty.tests.data_dir, 'vcf2.vcf')
  db_name = os.path.join(mitty.tests.data_dir, 't2.h5')
  open(vcf_name, 'w').write(_vcf)
  p = vcf2pop.vcf_to_pop(vcf_fname=vcf_name, pop_fname=db_name)

  vl = p.get_sample_variant_list_for_chromosome(1, 's1')
  assert vl[0][0]['pos'] == 119
  assert vl[0][0]['stop'] == 120
  assert vl[1][1]['pos'] == 139
  assert vl[1][1]['stop'] == 142

  assert vl[0].size == 2
  assert vl[1].size == 2

  db_name = os.path.join(mitty.tests.data_dir, 't3.h5')
  p = vcf2pop.vcf_to_pop(vcf_fname=vcf_name, pop_fname=db_name, sample_name='s2')

  vl = p.get_sample_variant_list_for_chromosome(1, 's2')
  assert vl[1][0]['pos'] == 119
  assert vl[1][0]['stop'] == 120

  assert vl[0].size == 1
  assert vl[1].size == 1