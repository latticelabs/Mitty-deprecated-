"""Convert a H5 file from one version to another"""
import click
import numpy as np
import h5py

str_dt = h5py.special_dtype(vlen=bytes)


def is_version_earlier_than(name, vs1):
  """Is vs0 < vs1 ?

  Expects vs0, vs1 in semantic versioning format X.Y.Z (possibly followed by .dev0)

  :param vs0:
  :param vs1:
  :return:
  """
  fp = h5py.File(name=name, mode='r')
  vs0 = fp.attrs['Mitty version']
  v0 = vs0.split('.')
  v1 = vs1.split('.')
  if v0[0] > v1[0]: return False  # 2.Y.Z vs 1.Y.Z
  if v0[0] < v1[0]: return True   # 1.Y.Z vs 2.Y.Z
  if v0[1] > v1[1]: return False  # 2.2.Z vs 2.1.Z
  if v0[1] < v1[1]: return True   # 2.1.Z vs 2.2.Z
  if v0[2] >= v1[2]: return False  # 2.2.Z vs 2.1.Z
  if v0[2] < v1[2]: return True   # 2.1.Z vs 2.2.Z


def migrate_pre1_34_0_to_1_34_0(db_name):
  """This modifies the database in place. You may want to make a back up."""
  fp = h5py.File(name=db_name, mode='r+')
  genome_metadata = []
  gmeta_keys = ['seq_id', 'seq_len', 'seq_md5']
  for k in [key for key in fp.keys() if key.startswith('chrom_')]:
    chrom = k[6:]
    print 'Processing chrom {}'.format(chrom)
    genome_metadata.append({gk: fp[k].attrs[gk] for gk in gmeta_keys})
    if 'master_list' in fp[k]:
      fp['/master_list/{}'.format(chrom)] = fp['{}/master_list'.format(k)]  # h5py.SoftLink('{}/master_list'.format(k))
    if 'samples' in fp[k]:
      for s in fp[k]['samples'].keys():
        fp['/samples/{}/{}'.format(s, chrom)] = fp['{}/samples/{}'.format(k, s)] # h5py.SoftLink('{}/samples/{}'.format(k, s))

  set_genome_metadata(fp, genome_metadata)
  fp.attrs['Mitty version'] = '1.34.0'


def set_genome_metadata(fp, genome_metadata):
  """Save chromosome sequence metadata

  :param genome_metadata: [{seq_id, seq_len, seq_md5} ...] in same order as seen in fa.gz file
                         same format as returned by Fasta.get_seq_metadata
  """
  dtype = [('seq_id', str_dt), ('seq_len', 'i4'), ('seq_md5', str_dt)]
  meta = [[gm[k] for gm in genome_metadata] for k in ['seq_id', 'seq_len', 'seq_md5']]
  fp.create_dataset('/ref_genome_meta', shape=(len(genome_metadata),), dtype=dtype,
                    data=np.core.records.fromarrays(meta, dtype))


@click.command()
@click.argument('name')
def cli(name):
  """This script modifies a pre 1.34.0 database to a 1.34.0 and later database"""
  if is_version_earlier_than(name, '1.34.0'):
    migrate_pre1_34_0_to_1_34_0(db_name=name)


if __name__ == '__main__':
  cli()