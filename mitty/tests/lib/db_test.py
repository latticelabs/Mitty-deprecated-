import tempfile
import os
from nose.tools import assert_list_equal, assert_dict_equal

from mitty.lib import db


def sample_roundtrip_test():
  """Sample roundtrip (save and load from database)"""
  temp_fp, temp_name = tempfile.mkstemp(suffix='.sqlite3')
  os.close(temp_fp)
  conn = db.connect(temp_name)

  c1 = [(1, 4, 'CAA', 'C', db.HOMOZYGOUS),
        (13, 14, 'C', 'G', db.HET_10),
        (20, 21, 'T', 'C', db.HET_10),
        (26, 27, 'T', 'TCGA', db.HOMOZYGOUS)]
  c2 = [(4, 6, 'CA', 'C', db.HOMOZYGOUS),
        (13, 14, 'C', 'G', db.HET_10),
        (20, 21, 'T', 'C', db.HET_10),
        (26, 27, 'T', 'TCGA', db.HOMOZYGOUS)]
  g = {'1': c1, '2': c2}

  db.save_sample('s1', g, conn)
  rt_g = {'1': [], '2': []}
  db.load_sample('s1', rt_g, conn)
  assert_dict_equal(g, rt_g, rt_g)
  os.remove(temp_name)


def sample_overwrite_test():
  """Overwrite sample"""
  temp_fp, temp_name = tempfile.mkstemp(suffix='.sqlite3')
  os.close(temp_fp)
  conn = db.connect(temp_name)
  c1 = [(1, 4, 'CAA', 'C', db.HOMOZYGOUS),
        (13, 14, 'C', 'G', db.HET_10)]
  c2 = [(4, 6, 'CA', 'C', db.HOMOZYGOUS),
        (20, 21, 'T', 'C', db.HET_10)]
  g = {'1': c1}
  db.save_sample('s1', g, conn)
  g = {'2': c2}
  db.save_sample('s1', g, conn)
  rt_g = {'2': []}
  db.load_sample('s1', rt_g, conn)
  assert_dict_equal(g, rt_g, rt_g)
  os.remove(temp_name)
