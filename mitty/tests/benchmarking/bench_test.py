import tempfile
import os
import json
import hashlib

from mitty.benchmarking.bench import *


def test_tuplelist2meta():
  """Benchmarking: Metadata round trip via json"""
  fi1 = {
    "tag": "f1",
    "sample": "1"
  }

  fi2 = {
    "tag": "f2",
    "sample": "2",
    "heights": [1, 2, 3]
  }

  metadata = OrderedDict([
    ("bench_run", "R1"),
    ("bench_name", "B1"),
    ("inputs", OrderedDict([
      ("inpZ", fi1),
      ("inpA", fi2)
    ])),
    ("tool", "gral-0.1.0")
  ])

  assert tuplelist2meta(meta2tuplelist(metadata)) == metadata

  _, t_file = tempfile.mkstemp(dir='./')
  json.dump(meta2tuplelist(metadata), open(t_file, 'w'))
  meta_from_file = json.load(open(t_file, 'r'))
  os.remove(t_file)
  assert tuplelist2meta(meta_from_file) == metadata


def test_meta2filename():
  """Benchmarking: Metadata to file name"""
  fi1 = {
    "tag": "f1",
    "sample": "1"
  }

  fi2 = {
    "tag": "f2",
    "sample": "2",
    "heights": [1, 2, 3]
  }

  metadata = OrderedDict([
    ("bench_run", "R1"),
    ("bench_name", "B1"),
    ("inputs", OrderedDict([
      ("inpZ", fi1),
      ("inpA", fi2)
    ])),
    ("tool", "gral-0.1.0")
  ])

  correct_fn = 'R1.B1.inpZ-f1.inpA-f2.gral-0.1.0'
  hu_fn = create_filename_prefix_from_metadata(metadata, use_hash=False)
  assert hu_fn == correct_fn, hu_fn

  hash_fn = create_filename_prefix_from_metadata(metadata)
  assert hash_fn == hashlib.md5(correct_fn).hexdigest()