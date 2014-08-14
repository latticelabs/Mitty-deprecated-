from setuptools import setup, find_packages

setup(
    name='mitty',
    version='1.0.0',
    packages=find_packages(exclude=['SbgWrappers', 'README-DATA', 'BigData']),
    install_requires=[
      'docopt >= 0.6.1',
      'pysam >= 0.8.0',
      'pyvcf',
      'numpy',
      'matplotlib',
      'h5py',
      'nose',
      'scipy',
      'bitarray'
    ]
)