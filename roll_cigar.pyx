## This function is the bottleneck, taking 1ms to run per call
import numpy as np
cimport numpy as np
#DTYPE = np.uint
#ctypedef np.uint_t DTYPE_t
def roll_cigar(this_read, np.ndarray[np.uint32_t, ndim=1] p_arr):
  """
  You can 'read along' with these tests from the Readme developers section

  Test a fully matching read
  >>> t_read = ['ATTG','~~~~', 0]; \
  p_arr = [1, 2, 3, 4, 5, 5, 5, 6, 9]; \
  roll_cigar(t_read, p_arr)
  (1, '4M')

  Test for read with insert
  >>> t_read = ['TTGT', '~~~~', 1]; \
  roll_cigar(t_read, p_arr)
  (2, '3M1I')

  Another test for read with insert
  >>> t_read = ['TGTT', '~~~~', 2]; \
  roll_cigar(t_read, p_arr)
  (3, '2M2I')

  Test for read with delete at end - should not show up in CIGAR
  >>> t_read = ['TTAC', '~~~~', 4]; \
  roll_cigar(t_read, p_arr)
  (5, '2I2M')

  Test for read spanning a deletion - should get a delete
  >>> t_read = ['ACAC', '~~~~', 0]; \
  p_arr = [1, 2, 5, 6, 7, 8, 9]; \
  roll_cigar(t_read, p_arr)
  (1, '2M2D2M')

  We actually missed this case: read with one matching base and then a delete
  >>> t_read = ['CACT', '~~~~', 1]; \
  roll_cigar(t_read, p_arr)
  (2, '1M2D3M')

  Test for an unmapped read: pos and cigars should be None
  >>> t_read = ['AATT', '~~~~', 2]; \
  p_arr = [1, 2, 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9]; \
  roll_cigar(t_read, p_arr)
  (0, '')

  """
  cdef int coord, counter, dp
  mapped = False
  cigar = ''
  coord = this_read[2]
  counter = 0
  cigar_fragment = None
  for n in range(coord, coord + len(this_read[0])):
    dp = p_arr[n+1] - p_arr[n]
    if dp == 1:
      mapped = True  # As long as we have one M we are a mapped read
      if cigar_fragment != 'M':
        if counter > 0:  # Flush
          cigar += '{:d}{:s}'.format(counter, cigar_fragment)
          counter = 0
      cigar_fragment = 'M'
      counter += 1
    elif dp == 0:
      if cigar_fragment != 'I':
        if counter > 0:  # Flush
          cigar += '{:d}{:s}'.format(counter, cigar_fragment)
          counter = 0
      cigar_fragment = 'I'
      counter += 1
    elif dp > 1:
      mapped = True  # As long as we have one M we are a mapped read
      if cigar_fragment != 'M':
        if counter > 0:  # Flush
          cigar += '{:d}{:s}'.format(counter, cigar_fragment)
          counter = 0
      cigar_fragment = 'M'  # We need to set this because we could be at the start of a read and type = None still
      counter += 1
      cigar += '{:d}{:s}'.format(counter, cigar_fragment)
      cigar_fragment = 'D'
      counter = dp - 1

  if cigar_fragment != 'D':  # Flush all but 'D'. We only write D if we cross a D boundary
    cigar += '{:d}{:s}'.format(counter, cigar_fragment)

  if mapped:
    align_pos = p_arr[coord]
  else:
    align_pos = 0
    cigar = ''

  return align_pos, cigar
