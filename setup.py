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
    install_requires=[
        'cython-ext',
        'numpy>=1.9.0',
        'scipy>=0.14.0',
        'docopt>=0.6.2',
        'pysam>=0.8.0',
        'PyVCF'],
    cython_ext='mitty/lib/*.pyx'
)
