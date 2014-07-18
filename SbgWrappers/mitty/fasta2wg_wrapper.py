"""This is the wrapper for fasta2wg. It also adds the feature of saving a concatenated .fa.gz file that can be used by
bwa etc.

Command line parameters are

fasta2wg  --index=IDX  --wg=WG  [-v]

"""
from sbgsdk import define, Process
import os
import json


#@require(mem_mb=4048, cpu=require.CPU_SINGLE)
class WholeGenomeCompactor(define.Wrapper):
  class Inputs(define.Inputs):
    fasta_list = define.input(name='FASTA files', description='Input .fa.gz files. One sequence per file.',
                              required=True, list=True)
    index = define.input(name='Index',
                         description='Index file (.json) that defines which file is which chromosome etc.',
                         required=True)

  class Outputs(define.Outputs):
    wg = define.output(name='.h5', description='Whole genome HDF5 file')
    fasta = define.output(name='fa.gz', description='Whole genome fasta file (can be sent to BWA etc.')

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

    with open(self.inputs.index, 'r') as f:
      index = json.load(f)
    # This is the only tricky thing: We need to modify the entered index file by replacing the bare name of the file with
    # the absolute path of the file.
    basename_list = [os.path.basename(fa) for fa in self.inputs.fasta_list]
    absolute_path_name_list = []
    for chrom_copy_list in index['chromosomes']:
      these_chrom_files = []
      for chrom_fname in chrom_copy_list:
        try:
          these_chrom_files.append(self.inputs.fasta_list[basename_list.index(chrom_fname)])
        except ValueError:
          print 'File in index missing from inputs'
      absolute_path_name_list.append(these_chrom_files)
    index['chromosomes'] = absolute_path_name_list

    index_file_name = os.path.join(output_dir, 'index_file.json')
    with open(index_file_name, 'w') as f:
      json.dump(index, f, indent=2)
    wg_file_name = os.path.join(output_dir, os.path.basename(os.path.splitext(self.inputs.index)[0]) + '.h5')
    fasta_file_name = os.path.join(output_dir, os.path.basename(os.path.splitext(self.inputs.index)[0]) + '.fa.gz')

    p = Process('python', '/Mitty/fasta2wg.py', '--index', index_file_name, '--wg', wg_file_name,
                '--fa', fasta_file_name, '-v')
    p.run()
    self.outputs.wg = wg_file_name
    self.outputs.wg.meta = self.outputs.wg.make_metadata(filetype='wg.gz')
    self.outputs.fasta = fasta_file_name
    self.outputs.fasta.meta = self.outputs.fasta.make_metadata(filetype='fa.gz')


def test_compactor():
  """."""
  json.dump(
    {
        "header": {
            "species": "Chimera Haploid Reference"
        },
        "chromosomes": [
            ["porcine_circovirus.fa.gz"],
            ["adenovirus.fa.gz"],
            ["herpes.fa.gz"],
            ["parvovirus.fa.gz"]
        ]
    }, open('/sbgenomics/test-data/wg_test.json', 'w'), indent=2)


  inputs = {
    'fasta_list': ['/sbgenomics/test-data/adenovirus.fa.gz',
                   '/sbgenomics/test-data/herpes.fa.gz',
                   '/sbgenomics/test-data/parvovirus.fa.gz',
                   '/sbgenomics/test-data/porcine_circovirus.fa.gz'],
    'index': '/sbgenomics/test-data/wg_test.json'
  }
  params = {}
  wrp = WholeGenomeCompactor(inputs, params)
  outputs = wrp.test()

  # Make sure the outputs exist
  assert os.path.exists(outputs.wg)
  assert os.path.exists(outputs.fasta)