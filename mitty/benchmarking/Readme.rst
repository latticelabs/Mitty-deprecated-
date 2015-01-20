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
1. Filesets are json files that aggregate collections of local files. This can be automated/regularized to be some 
kind of repository later to ensure versioning and uniqueness of files.
2. Tools are Python plugins that are wrappers around locally installed tools accessible via the command line. The wrapper 
constructs the command line required to run the tool and runs any appropriate pre- and post-processing to make the tool
run
3. Tool parameter sets are also json files listing a particular set of parameters for that tool
4. A benchmark parameter set is a json file that lists a particular set of parameters for a benchmark criterion
5. A benchmark file is another simple json file that lists the filesets, tools, tool parameter sets and bench mark parameter
sets that go into a benchmark. The benchmarker loads the benchmark spec .json file and runs the benchmarks

6. A to be determined visualizer then displays visualizations of the tool performance

Why not do this via the platform?
--------------------------------
The platform does not have loops or other constructs necessary to make this kind of run convenient. In order to run this
benchmark we would need to construct multiple pipelines, one for each parameter combination, which is tedious and error
prone.


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