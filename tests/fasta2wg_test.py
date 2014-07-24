import os
import tempfile
import h5py
from mitty.fasta2wg import save_genome_to_hdf5

fname = tempfile.gettempdir() + '/test.h5'


def write_test():
  """Testing creation of WG file from index and verifying its contents."""

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
  with h5py.File(fname, 'w') as h5_fp:
    save_genome_to_hdf5(index, h5_fp)

  #Now check to see that it exists
  assert os.path.exists(fname)

  #Now check some details of the file that convinces us fasta2wg is working as expected
  with h5py.File(fname, 'r') as h5_fp:
    assert h5_fp['sequence/2/1'].size == h5_fp['sequence/2/2'].size
    assert h5_fp['sequence/3/1'].attrs['reference']

