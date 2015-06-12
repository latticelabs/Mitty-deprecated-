from setuptools import setup, find_packages


setup(
    name='mitty',
    version='1.3.3.dev',
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
      # Register example tool wrapper
      'mitty.benchmarking.tools': ['bwa = mitty.benchmarking.tool_wrappers.bwa'],
      # Command line scripts
      'console_scripts': ['genomes = mitty.genomes:cli',
                          'reads = mitty.reads:cli',
                          'plot_align = mitty.util.plot_align:cli [mplot]',
                          'plot_gc_bias = mitty.util.plot_gc_bias:cli [mplot]',
                          'perfectbam = mitty.util.perfectbam:cli',
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
      'pysam>=0.8.1'
    ],
    extras_require = {
      'mplot': ['matplotlib>=1.3.0']
    }
)