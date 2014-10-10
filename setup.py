from setuptools import setup, find_packages

setup(
    name='mitty',
    version='1.0.0',
    description='Simulator for genomic data',
    author='Seven Bridges Genomics',
    author_email='kaushik.ghose@sbgenomics.com',
    packages=find_packages(include=['mitty*'], exclude=['mitty.util']),
    scripts=['mitty/denovo.py', 'mitty/checkbam.py', 'mitty/population.py', 'mitty/reads2bam.py', 'mitty/vcf2reads.py'],
    # install_requires=[
    #   'distribute>=0.7.3',
    #   'setuptools>=0.7',
    #   'Cython==0.21',
    #   'PyVCF>=0.6.7',
    #   'docopt==0.6.2',
    #   'mock==1.0.1',
    #   'nose==1.3.4',
    #   'numpy==1.9.0',
    #   'pyparsing==2.0.2',
    #   'pysam==0.8.0',
    #   'python-dateutil==2.2',
    #   'six==1.7.3',
    #   'wsgiref==0.1.2'
    # ],
    # dependency_links=['git+https://github.com/jamescasbon/PyVCF.git'],
    include_package_data=True,
    package_data={'mitty': ['tests/data/*']}
)
