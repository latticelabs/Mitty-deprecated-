import sys
sys.path.append('../mitty')  # This lets tests be run from pretty much anywhere

import tempfile
import os
import h5py
from mitty.fasta2wg import save_genome_to_hdf5
from shutil import rmtree

# These need to be available to the rest of the test suite for the mutations module
data_dir = tempfile.gettempdir()
wg_name = os.path.join(data_dir, 'test.h5')


def setup_package():
  """In order to speed up tests we create a complete chain of data starting from a whole genome file and ending at
  simulated reads. If this function fails it means Mitty is broken in some fundamental way."""
  if not os.path.exists(data_dir):
    os.makedirs(data_dir)

  # This setup only works when called from the project root directory such that Data is accessible
  index = {
    "header": {
        "species": "Test Chimera"
    },
    "chromosomes": [
        ["Data/adenovirus.fa.gz"],
        ["Data/porcine_circovirus.fa.gz", "Data/porcine_circovirus.fa.gz"],
        ["Data/altered_porcine.fa.gz"]
    ]
  }
  with h5py.File(wg_name, 'w') as h5_fp:
    save_genome_to_hdf5(index, h5_fp)


def teardown_package():
  rmtree(data_dir)