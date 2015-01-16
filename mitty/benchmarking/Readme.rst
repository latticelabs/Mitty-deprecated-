Goals
=====

Command line tool that will input a benchmark-suite file, run benchmarks declared in the file and then produce a result 
directory and a javascript-ed webpage that allows interactive exploration of the results.

Specs
=====
* Simulated ground truth data sets at different levels
  * HG38, SNPs, indels 
  * Perfect, corrupted
* Each benchmark consists of
  * A set of tools and a set of parameters
  * A set of input files

* A benchmark run iterates over all sets, skipping any combinations so marked
* Data files can be local or located in a git repository
* Results are placed in subdirectories based on combinations
* A d3 driven data explorer
  * Drill down of results
  


Architecture
============
1. Tools and file sets are Python plugins. file sets are packages containing data and a small Python class that 
describes them. ``pip install``-ing the packages results in them being registered as ``mitty.benchmarking`` components
which can then be looked up. This has the advantage that the data can be pulled from a central repository, can be 
versioned etc. A file set can depend on another fileset, such that we don't need to replicate common elements. For example
the hg38 can be a file set. This, by itself, is not useful, but another file set containing null reads from hg38 would
have it as a dependency (since aligners need the reference). Another file set containing a VCF referenced to hg38 and 
reads could also have it as a dependency, so that installing these sets does not duplicate hg38


2. A simple command-line application lets us list the tools and file sets available and add them to the benchmark. Then,
based on the tools selected it guides us through setting parameter sets and adding benchmark parameter sets. It then
saves a .json file with the benchmark specs. This .json file can be written by hand or any other tool too.

3. The benchmarker loads the benchmark spec .json file and runs the benchmarks

4. A to be determined visualizer then displays visualizations of the tool performance

Why not do this via the platform?
--------------------------------
The platform does not have loops or other constructs necessary to make this kind of run convenient. In order to run this
benchmark we would need to construct multiple pipelines, one for each parameter combination, which is tedious and error
prone.







---

1. A benchmark-suite file is 

1. A Python wrapper is written for the tool that acts as an interface between the benchmarker and the tool. The wrapper 
constructs the command line required to run the tool and runs any appropriate pre- and post-processing to make the tool
run
2. A parameter file tells the bench-marker which input file suites to use, what the parameter grid is and what wrapper
to use
3. The output is displayed as an interactive web-page that allows you to drill down the results.


Result organization
===================

+------------+-------------+-------------+ 
|            | Tool 1      | Tool 2      | 
+------------+------+------+------+------+
|            |  p1  |  p2  |  p1  |  p2  |
+============+======+======+======+======+
| Set 1      |      |   X  |      |  X   | 
+------------+------+------+------+------+ 
| Set 2      |      |      |  X   |      | 
+------------+------+------+------+------+ 
 
This grid is repeated for each bench marker parameter set (The benchmark itself can have tolerances etc. which don't
change the tool results but change the assessment of the results)


Result files are stored in separate directories

   |---------------- file set
   | |-------------- bench mark parameter
   | | |------------ tool 
   | | | |---------- tool parameter set
   | | | |
  s1bp1t1p1



Test Data
=========
The test data are stored on a git repository, say github, 




Interactive data explorer: alignments
=====================================
* Confusion matrix showing 2d histogram of where reads end up
* 
* Variations 

* Chromosomes are displayed vertically, 


TODO
====
Think about how to compare VCFs