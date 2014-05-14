from setuptools import setup, find_packages

setup(
    name='mitty',
    version='1.0.0',
    packages=find_packages(exclude=['SbgWrappers', 'test-data']),
    install_requires=['docopt, pysam, pyvcf, numpy']
)
