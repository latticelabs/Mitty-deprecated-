from setuptools import setup, find_packages

setup(
    name='mitty',
    version='1.0.0',
    description='Simulator for genomic data',
    author='Seven Bridges Genomics',
    author_email='kaushik.ghose@sbgenomics.com',
    packages=find_packages(include=['mitty*'], exclude=['mitty.util']),
    scripts=['mitty/denovo.py', 'mitty/checkbam.py', 'mitty/population.py', 'mitty/reads2bam.py', 'mitty/vcf2reads.py'],
    include_package_data=True,
    package_data={'mitty': ['tests/data/*']},
    setup_requires=['cython_ext'],
    install_requires=[
        'setuptools>=7.0',
        'numpy>=1.9.0',
        'scipy>=0.14.0',
        'docopt>=0.6.2',
        'pysam>=0.8.0',
        'PyVCF'],
    dependency_links=['git+https://gitlab.sbgenomics.com:9443/bjoern/cython_ext.git#egg=cython_ext-0.1',
                      'git+https://github.com/jamescasbon/PyVCF.git#egg=PyVCF'],
    cython_ext='mitty/lib/*.pyx'
)
