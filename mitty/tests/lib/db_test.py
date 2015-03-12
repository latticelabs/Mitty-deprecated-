import tempfile
import os
from nose.tools import assert_list_equal, assert_dict_equal

from mitty.lib import db
from mitty.lib.variation import HOMOZYGOUS as HOM, HET_01, HET_10


def sample_roundtrip_test():
  """Sample roundtrip (save and load from database)"""
  temp_fp, temp_name = tempfile.mkstemp(suffix='.sqlite3')
  os.close(temp_fp)
  conn = db.connect(temp_name)

  c1 = [(1, 4, HOM, 'CAA', 'C'),
        (13, 14, HET_10, 'C', 'G'),
        (20, 21, HET_10, 'T', 'C'),
        (26, 27, HOM, 'T', 'TCGA')]
  c2 = [(4, 6, db.HOMOZYGOUS, 'CA', 'C'),
        (13, 14, db.HET_10, 'C', 'G'),
        (20, 21, db.HET_10, 'T', 'C'),
        (26, 27, db.HOMOZYGOUS, 'T', 'TCGA')]
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
  c1 = [(1, 4, db.HOMOZYGOUS, 'CAA', 'C'),
        (13, 14, db.HET_10, 'C', 'G')]
  c2 = [(4, 6, db.HOMOZYGOUS, 'CA', 'C'),
        (20, 21, db.HET_10, 'T', 'C')]
  g = {'1': c1}
  db.save_sample('s1', g, conn)
  g = {'2': c2}
  db.save_sample('s1', g, conn)
  rt_g = {'2': []}
  db.load_sample('s1', rt_g, conn)
  assert_dict_equal(g, rt_g, rt_g)
  os.remove(temp_name)
