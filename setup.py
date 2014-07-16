from setuptools import setup, find_packages

setup(
    name='mitty',
    version='1.0.0',
    packages=find_packages(exclude=['SbgWrappers', 'README-DATA', 'BigData']),
    install_requires=['docopt', 'pysam', 'pyvcf', 'numpy', 'matplotlib']
)