from setuptools import setup

setup(
    name='mitty',
    version='1.0.0',
    description='Simulator for genomic data',
    author='Seven Bridges Genomics',
    author_email='kaushik.ghose@sbgenomics.com',
    packages=['mitty', 'mitty.lib', 'mitty.plugins', 'mitty.plugins.reads', 'mitty.plugins.variants', 'mitty.util'],
    scripts=['mitty/denovo.py', 'mitty/checkbam.py', 'mitty/population.py', 'mitty/reads2bam.py', 'mitty/vcf2reads.py'],
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