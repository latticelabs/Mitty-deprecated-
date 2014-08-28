import tempfile
import os
import h5py
from mitty.fasta2wg import save_genome_to_hdf5
import mitty.variation
from shutil import rmtree

# These need to be available to the rest of the test suite
source_tree_root = os.path.join(os.path.dirname(__file__), os.path.pardir)
example_fasta_genome = os.path.join(source_tree_root, 'examples', 'data')
data_dir = 'mitty_test_data_dir'  # tempfile.mkdtemp()
small_vcf_name = os.path.join(data_dir, 'small.vcf')


class SimulatedHDF:
  """We needed a class that pretended it is a HDF5 file for some of our tests."""

  def __init__(self):
    self.data = {}

  def insert(self, path, value):
    keys = path.split('/')
    this = self.data
    for key in keys[:-1]:
      if key == '': continue
      try:
        this = this[key]
      except KeyError:
        this[key] = {}
        this = this[key]
    this[keys[-1]] = value

  def __getitem__(self, item):
    keys = item.split('/')
    this = self.data
    for key in keys:
      if key == '': continue
      this = this[key]  # Custom error handling here?
    return this


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


def create_small_vcf():
  with open(small_vcf_name, 'w') as fp:
    fp.write(
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
      "1\t1\t.\tC\tCAA\t100\tPASS\t.\tGT\t0/1\n"
      "1\t3\t.\tCAG\tC\t100\tPASS\t.\tGT\t1/0\n"
      "2\t7\t.\tG\tT\t100\tPASS\t.\tGT\t0/1\n"
      "2\t9\t.\tGTT\t.\t100\tPASS\t.\tGT\t1/0\n"
      "2\t13\t.\tGTT\tTTG\t100\tPASS\t.\tGT\t1/1\n"
    )
  mitty.variation.compress_and_index_vcf(small_vcf_name, small_vcf_name + '.gz')


def reset_model(mod):
  vg = mod['model'].variant_generator
  for d in vg.__dict__.keys():
    delattr(vg, d)


def setup_package():
  """In order to speed up tests we create a complete chain of data starting from a whole genome file and ending at
  simulated reads. If this function fails it means Mitty is broken in some fundamental way."""
  os.makedirs(data_dir)
  #create_wg()
  create_small_vcf()


def teardown_package():
  rmtree(data_dir)
