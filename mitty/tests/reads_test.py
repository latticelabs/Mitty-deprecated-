import os
import json

from click.testing import CliRunner

import mitty.lib.mio as mio
import mitty.lib.variants as vr
import mitty.reads as reads
import mitty.tests


def integration_test():
  """'reads' command line program"""
  param_file = os.path.abspath(os.path.join(mitty.tests.data_dir, 'param.json'))
  db_file = os.path.abspath(os.path.join(mitty.tests.data_dir, 'pop.hdf5'))
  read_prefix = os.path.abspath(os.path.join(mitty.tests.data_dir, 'reads'))

  test_params = {
    "files": {
      "reference_dir": mitty.tests.example_data_dir,  # Use this if the reference consists of multiple .fa files in a directory
      "dbfile": db_file,
      "output_prefix": read_prefix,
      "interleaved": True
    },
    "sample_name": "g0_s0",   # Name of sample
    "rng": {
      "master_seed": 1
    },
    "chromosomes": [1, 2],               # List the chromosomes the reads should be taken from
    "variants_only": False,             # If true, reads will only come from the vicinity of variants
    "corrupt": True,                    # If true, corrupted reads will also be written
    "coverage": 1,                     # Coverage
    "coverage_per_block": 0.1,
    "read_model": "simple_illumina",    # Model specific parameters, need to be under the key "model_params"
    "model_params": {
      "read_len": 100,          # length of each read
      "template_len_mean": 250, # mean length of template
      "template_len_sd": 30,    # sd of template length
      "max_p_error": 0.01,      # Maximum error rate at tip of read
      "k": 20                   # exponent
    }
  }

  json.dump(test_params, open(param_file, 'w'))

  r_seq = mio.Fasta(multi_dir=mitty.tests.example_data_dir)

  pos = [27]
  stop = [28]
  ref = ['T']
  alt = ['G']
  p = [0.9]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  ml.sort()
  index_list = vr.l2ca([(0, 2)])
  pl = vr.Population(fname=db_file, mode='w', in_memory=False,
                     genome_metadata=r_seq.get_seq_metadata())
  pl.set_master_list(chrom=1, master_list=ml)
  pl.add_sample_chromosome(chrom=1, sample_name='g0_s0', indexes=index_list)

  runner = CliRunner()
  result = runner.invoke(reads.cli, ['generate', param_file])
  assert result.exit_code == 0, result
  assert os.path.exists(read_prefix + '.fq')
  assert os.path.exists(read_prefix + '_c.fq')