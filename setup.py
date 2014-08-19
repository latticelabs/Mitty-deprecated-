from setuptools import setup, find_packages

setup(
    name='mitty',
    version='1.0.0',
    packages=['mitty'],
    install_requires=[
      'docopt >= 0.6.1',
      'pysam >= 0.8.0',
      'pyvcf',
      'numpy',
      'matplotlib',
      'h5py',
      'nose',
      'scipy'
    ]
)