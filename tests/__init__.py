import tempfile
import os
import h5py
from mitty.fasta2wg import save_genome_to_hdf5
from shutil import rmtree

# These need to be available to the rest of the test suite for the mutations module
source_tree_root = os.path.join(os.path.dirname(__file__), os.path.pardir)
data_dir = tempfile.mkdtemp()
wg_name = os.path.join(data_dir, 'test.h5')


def create_wg():
  index = {
    "header": {
        "species": "Test Chimera"
    },
    "chromosomes": [
        [source_tree_root + "/Data/adenovirus.fa.gz"],
        [source_tree_root + "/Data/porcine_circovirus.fa.gz", source_tree_root + "/Data/porcine_circovirus.fa.gz"],
        [source_tree_root + "/Data/altered_porcine.fa.gz".format(source_tree_root)]
    ]
  }
  with h5py.File(wg_name, 'w') as h5_fp:
    save_genome_to_hdf5(index, h5_fp)


def setup_package():
  """In order to speed up tests we create a complete chain of data starting from a whole genome file and ending at
  simulated reads. If this function fails it means Mitty is broken in some fundamental way."""
  if not os.path.exists(data_dir):
    os.makedirs(data_dir)
  create_wg()


def teardown_package():
  rmtree(data_dir)
