.. _additional-utilities:

..
    Additional Utilities
    ====================

    `plot_gc_bias`: Plot GC bias of a BAM file
    ------------------------------------------
    There is a utility for inspecting the GC bias in BAM files called `plot_gc_bias`.

    .. command-output::  ../examples/reads_gc_bias/doc_script.sh 1> /dev/null  # A script to generate a BAM file with a GC bias
      :shell:



    (You will, of course, have noticed that the stock `simple_illumina` plugin has a GC bias parameter.)

    .. command-output:: plot_gc_bias ../examples/data/red_alga.fa.gz  ../examples/reads_gc_bias/Out/reads.bam --win 1000 --g0 0.45 --g1 0.65 --out ../examples/reads_gc_bias/Out/gc_bias.png
    .. image:: ../examples/reads_gc_bias/Out/gc_bias.png
