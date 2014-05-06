Examples
========

Complete example
-------------

We will start with `complete_example.sh` which is a short bash script that takes Mitty through its paces. This example
uses the porcine_circovirus sequence (included in the data directory) which is 702 bases long.



Note, how in `read_par.json` we don't have any parameters for read corruption. In such a case Mitty uses the default
parameters set in the code. If there are none, an error will result. In all the other examples we explicitly pass
all parameters as required as that is good practice.

