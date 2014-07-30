import json
import vcf
from nose.tools import assert_equals, ok_
from . import *
from mitty import mutate


def test_script():
  """mutate command line program"""

  ok_(os.path.exists(wg_name),
      msg='No whole genome file ({:s}). This should be created by package test setup in tests/__init__.py'.format(wg_name))

  tempdir = tempfile.mkdtemp()
  vcf_name = os.path.join(tempdir, 'mutate_test.vcf.gz')
  param_name = os.path.join(tempdir, 'mutate_test_param.json')

  args = {
    '--wg': wg_name,  # Our package wide setup as generated this
    '--vcf': vcf_name,
    '--paramfile': param_name,
    '--master_seed': None
  }

  # Test with a simple stock plugin
  json.dump({
    "variant_models": [
      {
        "chromosome": [2],  # Shortest chromosome in chimera
        "model": "snp",
        "phet": 0.5,
        "p": 0.01,
        "poisson_rng_seed": 1,
        "base_sub_rng_seed": 2
      }
    ]
  }, open(param_name, 'w'))

  mutate.main(args)

  """
  The output file should be (with variations in the header date and possibly command line)
  ##fileformat=VCFv4.1
  ##fileDate= ...
  ##source= ...
  ##reference=test.h5
  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample
  2       3       .       G       A       100     PASS    .       GT      1/1
  2       99      .       C       A       100     PASS    .       GT      1/1
  2       251     .       G       A       100     PASS    .       GT      0/1
  2       444     .       T       A       100     PASS    .       GT      1/1
  2       488     .       C       T       100     PASS    .       GT      1/1
  2       570     .       T       C       100     PASS    .       GT      1/1
  2       678     .       T       A       100     PASS    .       GT      0/1
  2       691     .       C       A       100     PASS    .       GT      0/1
  """
  ok_(os.path.exists(vcf_name), msg='Output file was not created')
  variants = [var for var in vcf.Reader(filename=vcf_name)]
  ok_(len(variants) == 8, msg='In correct number of variants')
  ok_(variants[7].REF == 'C', msg='File contents not as expected')

  rmtree(tempdir)  # Be neat