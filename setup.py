from setuptools import setup, find_packages

setup(
    name='mitty',
    version='1.0.0',
    description='Simulator for genomic data',
    author='Seven Bridges Genomics',
    author_email='kaushik.ghose@sbgenomics.com',
    packages=find_packages(include=['mitty*'], exclude=['mitty.util']),
    include_package_data=True,
    package_data={'mitty': ['tests/data/*']},
    entry_points={
      # Register the built in plugins
      'mitty.plugins.variants': ['snp = mitty.plugins.variants.snp_plugin'],
      'mitty.plugins.reads': ['simple_sequential = mitty.plugins.reads.simple_sequential_plugin',
                              'simple_illumina = mitty.plugins.reads.simple_illumina_plugin'],
      # Command line scripts
      'console_scripts': ['denovo = mitty.denovo:cli', 'vcf2reads = mitty.vcf2reads:cli']
    },
    install_requires=[
      'cython-ext',
      'numpy>=1.9.0',
      'scipy>=0.14.0',
      'docopt>=0.6.2',
      'pysam>=0.8.0',
      'PyVCF'
    ],
    cython_ext='mitty/lib/*.pyx'
)