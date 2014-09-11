from setuptools import setup, find_packages

setup(
    name='mitty',
    version='1.0.0',
    description='Simulator for genomic data',
    author='Seven Bridges Genomics',
    author_email='kaushik.ghose@sbgenomics.com',
    #packages=['mitty', 'mitty.lib', 'mitty.plugins', 'mitty.plugins.reads', 'mitty.plugins.variants', 'mitty.util'],
#    packages=['mitty', 'mitty.lib', 'mitty.plugins', 'mitty.plugins.reads', 'mitty.plugins.variants',
#              'tests', 'tests.'],
    packages=find_packages(include=['mitty', 'tests']),
    scripts=['mitty/denovo.py', 'mitty/checkbam.py', 'mitty/population.py', 'mitty/reads2bam.py', 'mitty/vcf2reads.py'],
    install_requires=[
      'distribute>=0.7.3',
      'setuptools>=0.7',
      'Cython==0.21',
      'PyVCF>=0.6.7',
      'docopt==0.6.2',
      'matplotlib==1.4.0',
      'mock==1.0.1',
      'nose==1.3.4',
      'numpy==1.9.0',
      'pyparsing==2.0.2',
      'pysam==0.8.0',
      'python-dateutil==2.2',
      'six==1.7.3',
      'wsgiref==0.1.2'
    ],
    data_files=[('examples/data', ['examples/data/chr1.fa',
                                   'examples/data/chr2.fa',
                                   'examples/data/chr3.fa',
                                   'examples/data/chr4.fa',
                                   'examples/data/chimera.fa',
                                   'examples/data/chimera.fa.gz'])]
)