import mitty.lib.vcf2db as vcf2db
import mitty.lib.db as mdb
import mitty.lib.variants as vr
import mitty.genomes as genomes

import os
import tempfile

seq_meta_data = [
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
  master_lists = [vr.VariantList(vd['pos'], vd['stop'], vd['ref'], vd['alt'], vd['p']) for vd in variant_data]

  conn1 = mdb.connect(db_name=':memory:')
  for n, meta in enumerate(seq_meta_data):
    mdb.save_chromosome_metadata(conn1, n + 1, **meta)

  for n, ml in enumerate(master_lists):
    ml.sort()
    mdb.save_master_list(conn1, n + 1, ml)

  for n, gt in enumerate(genotype_data):
    mdb.save_sample(conn1, 0, 0, n + 1, gt)

  temp_vcf_prefix = os.path.join(tempfile.gettempdir(), 'test')
  genomes.write_sample_vcfs(conn1, [0], temp_vcf_prefix)

  _, temp_db_file = tempfile.mkstemp(suffix='db')

  vcf2db.vcf_to_db(temp_vcf_prefix + '_s0.vcf', temp_db_file)
  conn2 = mdb.connect(db_name=temp_db_file)

  loaded_mls = [mdb.load_master_list(conn2, n) for n in [1, 2]]
  loaded_chroms = [mdb.load_sample(conn2, 0, 0, n) for n in [1, 2]]

  for n in [0, 1]:  # Chromosomes
    for v1, v2 in zip(loaded_chroms[n], genotype_data[n]):
      assert v1[1] == v2[1]  # Genotypes match
      for k in ['pos', 'stop', 'ref', 'alt']:
        assert loaded_mls[n].variants[v1[0]][k] == master_lists[n].variants[v2[0]][k]  # Variant data match

  os.remove(temp_vcf_prefix + '_s0.vcf')
  os.remove(temp_db_file)