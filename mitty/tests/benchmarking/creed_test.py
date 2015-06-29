from nose.tools import assert_raises

import mitty.benchmarking.creed as creed


def categorized_reads_class_test():
  """benchmarking: CategorizedReads class"""
  assert_raises(RuntimeError, creed.CategorizedReads)  # Need seq len and seq id for new file
  seq_names = ['chr1', 'chr2']
  seq_lengths = [100, 200]
  cat_reads = creed.CategorizedReads(seq_names=seq_names, seq_lengths=seq_lengths, copies=2)
  assert cat_reads.get_sequence_info()[1][1] == seq_lengths[1]

  cat_reads.append(1, 0, 20, 100, 0b10010)
  cat_reads.append(1, 0, 10, 100, 0b10000)
  cat_reads.append(0, 1, 10, 100, 0b10000)
  cat_reads.finalize()

  dta = cat_reads.get_data(1, 0)
  assert dta.shape[0] == 2
  assert dta['pos'][0] < dta['pos'][1]

  assert_raises(RuntimeError, cat_reads.append, 1, 0, 10, 100, 0b10000)  # Can't append after finalizing
