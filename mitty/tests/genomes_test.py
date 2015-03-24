from nose.tools import assert_raises

import os
import tempfile

import sqlite3 as sq

import mitty.lib.io as mio
import mitty.lib.db as mdb
import mitty.genomes as genomes
import mitty.plugins.variants.snp_plugin as snp
import mitty.plugins.site_frequency.double_exp as double_exp
from mitty.tests import *  # To get definitions from the setup script


def run_simulations_test():
  """Core simulation loop"""
  temp_fp, temp_name = tempfile.mkstemp(suffix='.sqlite3')
  os.close(temp_fp)

  pop_db_name = temp_name
  ref = mio.Fasta(multi_dir=example_data_dir)
  sfs_model = double_exp.Model()
  variant_models = [snp.Model(p=0.001), snp.Model(p=0.005)]
  chromosomes = [1, 2]
  sample_size = 20
  master_seed = 2

  genomes.run_simulations(pop_db_name, ref, sfs_model, variant_models, chromosomes, sample_size, master_seed)

  # The parts are tested separately, we just want a superficial test to see if we wrote the correct number of things
  conn = mdb.connect(temp_name)
  ml = mdb.load_master_list(conn, 1)
  assert len(ml) > 0
  ml = mdb.load_master_list(conn, 1)
  assert len(ml) > 0
  for n in range(sample_size):
    chrom = mdb.load_sample(conn, 0, n, 1)
    assert len(chrom) > 0
  for n in range(sample_size):
    chrom = mdb.load_sample(conn, 0, n, 2)
    assert len(chrom) > 0

  assert_raises(sq.OperationalError, mdb.load_sample, conn, 0, 0, 3)  # Should be no such table