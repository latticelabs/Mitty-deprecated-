from setuptools import setup, find_packages

__version__ = eval(open('mitty/version.py').read().split('=')[1])
setup(
    name='mitty',
    version=__version__,
    description='Simulator for genomic data',
    author='Seven Bridges Genomics',
    author_email='kaushik.ghose@sbgenomics.com',
    packages=find_packages(include=['mitty*']),
    include_package_data=True,
    entry_points={
      # Register the built in plugins
      'mitty.plugins.sfs': ['double_exp = mitty.plugins.site_frequency.double_exp'],
      'mitty.plugins.variants': ['snp = mitty.plugins.variants.snp_plugin',
                                 'delete = mitty.plugins.variants.delete_plugin',
                                 'uniformdel = mitty.plugins.variants.uniform_deletions',
                                 'uniformins = mitty.plugins.variants.uniform_insertions',
                                 'insert = mitty.plugins.variants.insert_plugin',
                                 #'inversion = mitty.plugins.variants.inversion_plugin',
                                 #'low_entropy_insert = mitty.plugins.variants.low_entropy_insert_plugin'
                                 ],
      'mitty.plugins.population': ['standard = mitty.plugins.population.standard'],
      'mitty.plugins.reads': ['simple_sequential = mitty.plugins.reads.simple_sequential_plugin',
                              'simple_illumina = mitty.plugins.reads.simple_illumina_plugin'],
      # Command line scripts
      'console_scripts': ['genomes = mitty.genomes:cli',
                          'reads = mitty.reads:cli',
                          'perfectbam = mitty.benchmarking.perfectbam:cli',
                          'alindel = mitty.benchmarking.indel_alignment_accuracy:cli',
                          'vcf2pop = mitty.lib.vcf2pop:cli',
                          'alindel_plot = mitty.benchmarking.indel_plot:cli',
                          'plot_align = mitty.util.plot_align:cli [plot]',
                          'plot_gc_bias = mitty.util.plot_gc_bias:cli [plot]',
                          'mismat = mitty.util.mismat:cli',
                          'splitta = mitty.util.splitta:cli',
                          'kmers = mitty.util.kmers:cli',
                          'pybwa = mitty.util.pybwa:cli']
    },
    install_requires=[
      'cython',
      'setuptools>=11.0.0',
      'numpy>=1.9.0',
      'docopt>=0.6.2',
      'click>=3.3',
      'pysam>=0.8.1',
      'h5py>=2.5.0'
    ],
    extras_require = {
      'plot': ['matplotlib>=1.3.0', 'scipy']
    }
)