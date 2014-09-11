from setuptools import setup, find_packages

setup(
    name='mitty',
    version='1.0.0',
    description='Simulator for genomic data',
    author='Seven Bridges Genomics',
    author_email='kaushik.ghose@sbgenomics.com',
    packages=['mitty', 'examples'],
    scripts=['denovo.py', 'checkbam.py', 'population.py', 'reads2bam.py', 'vcf2reads.py'],
    install_requires=[
      'docopt >= 0.6.1',
      'pysam >= 0.8.0',
      'pyvcf >= 0.6.7',
      'numpy >= 1.9.0',
      'matplotlib >= 1.4.0',
      'nose >= 1.3.4',
      'scipy >= 0.14.0'
    ]
)