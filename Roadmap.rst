v 3.0.0 (April 2016)

Core engine:
  * Computation by segments of sequence (rather than just by chromosome)
  * Parallelize most computations

Benchmarking: Core
  * Integrate benchmarking stuff with CWL
  * Fuzzer loop
  * Debug loop
    - trace wrong variant call back to reads on BAM


v 2.0.0 (Jan 2016)

Core engine:
  * Haplotype blocks
  * Repetitive regions can be ignored for taking reads
  * More efficient FASTA reading
  * Cythonize bottle necks

Benchmarking: Core
  * ... something to do with read qnames and encoding of the correct answer
  * Alignment accuracy for reads under sample variations
  * Alignment accuracy for reads under reference variations

Benchmarking: Visualizations
  * mis-alignment plots (circle and matrix)
  * vcf concordance plots
  * indel accuracy plots for alignment
  * indel accuracy plots for variant calling (TP)
  * ROC curves?

Documentation
  * Finish tutorial
    - gc bias section
  * Write up instructions on how to write plugins
  * Write up theory parts
    - site frequency spectrum, population models