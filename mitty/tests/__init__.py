import os
#import mitty.lib.variation
import mitty.lib.io
import mitty.lib
from shutil import rmtree
# These need to be available to the rest of the test suite
source_tree_root = os.path.join(os.path.dirname(__file__))
example_data_dir = os.path.join(source_tree_root, 'data')

example_fasta_genome = os.path.join(source_tree_root, 'data')
data_dir = 'mitty_test_data_dir'  # tempfile.mkdtemp()
small_vcf_name = os.path.join(data_dir, 'small.vcf')
null_fastq_name = os.path.join(data_dir, 'null_reads.fq')
fake_bam_name = os.path.join(data_dir, 'fake.bam')


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
  mitty.lib.io.compress_and_index_vcf(small_vcf_name, small_vcf_name + '.gz')


def create_null_reads():
  """Generate a file of paired reads from the chimera test genome."""
  import mitty.vcf2reads as vcf2reads
  from mitty.lib.genome import FastaGenome
  model_params = {
    "paired": True,
    "read_len": 100,
    "template_len": 250,
    "read_advance": 10
  }
  with open(null_fastq_name, 'w') as fastq_fp:
    ref = FastaGenome(example_fasta_genome, persist=True)
    read_model = mitty.lib.load_reads_plugin('simple_sequential')
    vcf2reads.main(fastq_fp, fastq_c_fp=None, ref=ref, g1={}, chrom_list=[3, 4],
                   read_model=read_model, model_params=model_params, block_len=10e6, master_seed=42)


def create_fake_bam():
  """Create a bam file with one correctly mapped read, one unmapped read and one incorrectly mapped read."""
  #mitty.lib.bam.sort_and_index_bam(fake_bam_name)


def setup_package():
  """In order to speed up tests we create a complete chain of data starting from a whole genome file and ending at
  simulated reads. If this function fails it means Mitty is broken in some fundamental way."""
  pass
  # os.makedirs(data_dir)
  # create_small_vcf()
  # create_null_reads()


def teardown_package():
  pass
  # rmtree(data_dir)