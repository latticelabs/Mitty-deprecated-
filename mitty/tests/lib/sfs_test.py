import mitty.lib.sfs as sitefs


def balancing_test():
  """Test balancing of site freq spectrum"""
  sfs = sitefs.QuantizedSfs(p=[float(n)/10.0 for n in range(10)], ideal_sfs=[0.1] * 10)
  iv = [[n]*n*1000 for n in range(10)]
  sfs.add(variants=iv)
  sfs.balance_sfs()
  assert len(sfs.quantized_master_list[0]) == 4500
  assert len(sfs.quantized_master_list[9]) == 4500



