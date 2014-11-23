from setuptools import setup, find_packages
from Cython.Build import cythonize

#extensions = [Extension("*", "*.pyx")]
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
    ext_modules=cythonize('mitty/lib/*.pyx'),
    entry_points={
      # Register the built in plugins
      'mitty.plugins.variants': ['snp = mitty.plugins.variants.snp_plugin'],
      # Command line scripts
      'console_scripts': ['denovo = mitty.denovo:cli']
    }
)