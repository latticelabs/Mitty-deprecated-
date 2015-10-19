from setuptools import setup, find_packages

setup(
    name='SimulationTools',
    version='1.19.0dev',
    packages=find_packages(exclude=['mock_data']),
    install_requires=['sbgsdk']
)