from setuptools import setup, find_packages
from Cython.Build import cythonize

setup(
    name='mitty',
    version='1.0.0',
    description='Simulator for genomic data',
    author='Seven Bridges Genomics',
    author_email='kaushik.ghose@sbgenomics.com',
    packages=find_packages(include=['mitty*']),
    include_package_data=True,
    entry_points={
      # Register the built in plugins
      'mitty.plugins.variants': ['snp = mitty.plugins.variants.snp_plugin',
                                 'delete = mitty.plugins.variants.delete_plugin',
                                 'bounded_delete = mitty.plugins.variants.bounded_len_delete_plugin',
                                 'insert = mitty.plugins.variants.insert_plugin',
                                 'inversion = mitty.plugins.variants.inversion_plugin',
                                 'low_entropy_insert = mitty.plugins.variants.low_entropy_insert_plugin'],
      'mitty.plugins.reads': ['simple_sequential = mitty.plugins.reads.simple_sequential_plugin',
                              'simple_illumina = mitty.plugins.reads.simple_illumina_plugin'],
      # Command line scripts
      'console_scripts': ['denovo = mitty.denovo:cli',
                          'vcf2reads = mitty.vcf2reads:cli',
                          'reads2bam = mitty.reads2bam:cli',
                          'splitta = mitty.util.splitta:cli',
                          'checkbam = mitty.benchmarking.checkbam:cli']
    },
    install_requires=[
      'setuptools>=0.7',
      'cython>=0.21.0',
      'numpy>=1.9.0',
      'scipy>=0.14.0',
      'docopt>=0.6.2',
      'pysam>=0.8.1',
      'PyVCF==0.7.0dev'
    ],
    ext_modules=cythonize('mitty/lib/*.pyx'),
)
