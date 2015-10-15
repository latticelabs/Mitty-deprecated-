"""For our tutorial we want to do demonstrations on a genome that
  * is not too big - don't want a ton of data to process each time we regenerate the docs or do a demo
  * want several chromosomes - important to illustrate all features
  * is not too small - we want a decent number of variations and reads

To this end we take the last five chromosomes of the red_alga genome and pretend these are chromosomes from a diploid
creature 'Reddus pentalgus'"""
import gzip

input_fname = 'red_alga.fa.gz'
output_fname = 'reddus_pentalgus.fa.gz'

skip_cntr = 16
with gzip.open(input_fname, 'r') as in_file, gzip.open(output_fname, 'w') as out_file:
  while skip_cntr > 0:
    ln = in_file.readline()
    if ln.startswith('>'):
      skip_cntr -= 1
  out_file.write(ln)
  for line in in_file:
    out_file.write(line)
