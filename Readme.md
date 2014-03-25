Mitty is a collection of modules and scripts that enable us to generate simulated genomic data to test our algorithms.


Dev notes
---------
### Parameter files for mutate
1. I chose to use parameter files because we often want to rerun experiments and it became clear early on that there would
be a lot of parameters.
1. I chose to use python for the parameter file for parsimony and flexibility
1. The parameter distribution between file and command line was based on predictions of which parameters we could
experiment with most during testing

(Since you were dying to know: Mitty comes from James Thurber's "The Secret Life of Walter Mitty" one of my favourite
pieces from one of my favourite authors. Though [Wikipedia][wiki] has a less favourable interpretation of what Walter Mitty
stands for I follow the interpretation found in the 2013 movie of the same name. Life is difficult and full of
insurmountable obstacles. If you do not even dream that you have surmounted these obstacles how are you going to even
start?)

[wiki]: http://en.wikipedia.org/wiki/The_Secret_Life_of_Walter_Mitty