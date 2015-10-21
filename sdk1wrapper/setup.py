from setuptools import setup, find_packages

setup(
    name='sbg_mitty',
    version='1.24.0dev0',
    packages=find_packages(exclude=['test-data']),
    install_requires=['sbgsdk']
)