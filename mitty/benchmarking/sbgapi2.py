"""
To use this set the CGC_TOKEN environment variable
"""

import os
import json
import time
from collections import OrderedDict
import logging

import requests

import bench

logger = logging.getLogger(__name__)

__base_url__ = 'https://cgc-api.sbgenomics.com/v2/'
__sbg_auth_token_env_var__ = 'CGC_TOKEN'


def join_url(*pieces):
  return '/'.join(s.strip('/') for s in pieces)


def get_auth_token(auth_token=None):
  """Convenient wrapper that either uses the passed auth token or loads it from an environment variable
  :param auth_token: authentication token to use
  """
  _auth_token = auth_token or os.environ.get(__sbg_auth_token_env_var__, None)
  if _auth_token is None:
    raise RuntimeError('SBG authentication token not found. Have you filled out {:s} ?'.format(__sbg_auth_token_env_var__))
  return _auth_token


# https://docs.sbgenomics.com/display/developerhub/Quickstart
def api(path='/', method='GET', query=None, data=None, auth_token=None):
  """Execute a query on the platform and return the response.

  :param path:
  :param method:
  :param query:
  :param data:
  :param auth_token:
  :return:
  """

  data = json.dumps(data) if isinstance(data, dict) else None
  base_url = __base_url__

  headers = {
    'X-SBG-Auth-Token': get_auth_token(auth_token),
    'Accept': 'application/json',
    'Content-type': 'application/json',
  }

  t0 = time.time()
  response = requests.request(method, join_url(base_url, path), params=query, data=data, headers=headers)
  response_dict = json.loads(response.content) if response.content else {}
  t1 = time.time()
  logger.debug('Request took {:f}s'.format(t1 - t0))

  if response.status_code / 100 != 2:
    print(response_dict)
    raise Exception('Server responded with status code %s.' % response.status_code)

  return response_dict


def get_billing_group(auth_token, billing_group_name=None):
  """Returns the first billing group if no name is supplied or name does not match"""
  billing_groups = api('/billing/groups', auth_token=auth_token)
  send_bills_to = billing_groups['items'][0]
  for grp in billing_groups['items']:
    if grp['name'] == billing_group_name:
      send_bills_to = grp
      break
  return send_bills_to


def get_existing_project(project_name, auth_token=None):
  p_list = [p for p in api('projects', auth_token=auth_token)['items'] if p['name'] == project_name]
  return p_list[0] if len(p_list) > 0 else None


def create_new_project(name, description=None, auth_token=None, billing_group_name=None):
  """Create new project. Return the details of the project

  :param name:
  :param description:
  :param auth_token:
  :param billing_group_name:
  :return:
  """
  _auth_token = get_auth_token(auth_token)
  billing_group = get_billing_group(auth_token=_auth_token, billing_group_name=billing_group_name)
  data = {"name": name, "description": description, "billing_group_id": billing_group['id']}
  return api('projects', method='POST', data=data, auth_token=_auth_token)


def get_or_create_new_project(name, description=None, auth_token=None, billing_group_name=None):
  """Create new project, or skip if project exists. Return the details of the project

  :param name:
  :param description:
  :param auth_token:
  :param billing_group_name:
  :return:
  """
  _auth_token = get_auth_token(auth_token)
  project = get_existing_project(name, _auth_token)
  return project or create_new_project(name, description, _auth_token, billing_group_name)


def delete_project(project, auth_token=None):
    return api('projects/{id:s}'.format(**project), method='DELETE', auth_token=auth_token)


# def get_credentials(project_name, user_name, auth_token=None):
#   """Setup a dictionary that carries info about our platform credentials, project etc."""
#   project_id = [project['id'] for project in api('projects', auth_token=auth_token)['items'] if project['name'] == project_name]
#   if len(project_id) == 0:
#     logger.error('No project named {:s}'.format(project_name))
#     return None
#   else:
#     project_id = project_id[0]
#   return {'AUTH_TOKEN': get_auth_token(auth_token), 'project_id': project_id, 'user': user_name}


def get_all_files_in_project(project, auth_token=None):
  """Returns the rich file info obtained from the platform, including id and metadata."""
  return api('projects/{id:s}/files'.format(**project), method='GET', auth_token=auth_token).get('items', [])


def get_file_by_name(f_info_l, file_name):
  for f in f_info_l:
    if f['name'] == file_name:
      return f
  else:
    return None


def copy_files(src_project, src_files, dest_project):
  f_list = get_all_files_in_project(src_project)
  return api('files/copy', method='POST', data={'project': dest_project['id'],
                                                'files': [f['id'] for f in f_list if f['name'] in src_files]})


def modify_file(project, original_file_name, new_file_name, new_metadata, src_f_list, auth_token=None):
  """Rename the file, setting the name and the metadata as required. Metadata is overwritten.
  No check is made to verify that the source f list is from the src project. That is the caller's responsibility

  :param project:
  :param original_file_name:
  :param new_file_name:
  :param new_metadata:
  :param src_f_list: If supplied, this saves us an API call to fetch the file list from the src_proj
  :param auth_token
  :return: The new file name details
  """
  f = get_file_by_name(src_f_list, original_file_name)
  if f is None:
    raise RuntimeError('File {:s} not found in file list'.format(original_file_name))

  request_body = {
    'name': new_file_name,
    'metadata': new_metadata
  }
  return api('files/{:s}'.format(f['id']), method='POST', data=request_body, auth_token=auth_token)


def get_all_apps_in_project(project, auth_token=None):
  """

  :param app_list:
  :param project:
  :param auth_token:
  :return:
  """
  return api('/apps/{id:s}'.format(**project)).get('items', [])  # Get all apps in the project


class SBGSDK2Executor(bench.BaseExecutor):
  """Execute benchmarking tasks on an SDK2 enabled platform"""

  def __init__(self, bench_run, project_name, auth_token=None):
    """Setup communication with the platform.

    :param bench_run:
    :param project_name:
    :param auth_token: Platform authentication token
    :return:
    """
    self.auth_token, self.bench_run, self.project = auth_token, bench_run, get_existing_project(project_name, auth_token)
    self._check_project()

    self.apps = get_all_apps_in_project(self.project, self.auth_token)
    self.files_in_project = get_all_files_in_project(self.project, self.auth_token)

    # Perform a few sanity checks
    run_name_ok, files_ok, apps_ok = self._check_run_name(), self._check_initial_files(), self._check_apps()
    # if not (run_name_ok and files_ok and apps_ok):
    #   raise RuntimeError('Error setting up benchmark. Please see log file')

  def _check_project(self):
    """Make sure project exists"""
    if self.project is None:
      raise RuntimeError('No project named {name:s} found on platform'.format(**self.project))

  def _check_run_name(self):
    """Make sure bench run name does not clash"""
    metadata = OrderedDict([
      ("bench_run", self.bench_run['bench_run_name']),
      ("bench_name", self.bench_run['bench_name']),
      ("bench_inputs", {})
    ])
    prefix = bench.create_filename_prefix_from_metadata(metadata, use_hash=self.bench_run['use_hash_for_filenames'])
    if len(filter(lambda x: x['name'].startswith(prefix), self.files_in_project)):
      logger.error('A bench run with this name already exists in the project. Please use a different name')
      return False
    return True

  def _check_initial_files(self):
    _file_names = [x['name'] for x in self.files_in_project]
    missing_files = filter(lambda x: x['file_name'] not in _file_names, self.bench_run['file_list'].itervalues())
    if len(missing_files):
      for f in missing_files:
        logger.error('File {:s} missing on platform project'.format(f['file_name']))
      return False
    return True

  def _check_apps(self):
    _app_names = [x['name'] for x in self.apps]
    missing_apps = filter(lambda x: x['app_name'] not in _app_names,
                          self.bench_run['tool_descriptions'].values() + self.bench_run['benchmark_tools'].values())
    if len(missing_apps):
      for f in missing_apps:
        logger.error('App {:s} missing on platform project'.format(f['app_name']))
      return False
    return True

  def start_job(self, app, input_files, output_files):
    """Find the app on the platform and start a new task with the given input files.

    :param app:
    :param input_files:
    :param output_files:
    :return:
    """
    # Refresh the file list
    self.files_in_project = get_all_files_in_project(self.project, self.auth_token)
    if not self._check_job_input_files(input_files):
      raise RuntimeError('Files required for {app_name:s} are missing'.format(**app))

    data = {
      "description": "Test",
      "name": "testsenad",
      "app_id": "Rfranklin/my-project/new-app/2",
      "project": "RFranklin/my-project",
      "inputs": {
        "my-input-node": {
          "class": "File",
          "path": "562785e6e4b00a1d67a8b1aa",
          "name": "example_human_known_indels.vcf"
        }
      }
    }

    app = find_app_in_project

    p = Process(target=self.test_process, args=(app, random.uniform(self.sleep_min, self.sleep_max), output_files))
    p.start()
    pid = str(p.pid)
    self.job_id[pid] = p
    self.job_files[pid] = {
      'app': app,
      'input_files': input_files,
      'output_files': output_files
    }
    return pid

  def _check_job_input_files(self, input_files):
    _file_names = [x['name'] for x in self.files_in_project]
    missing_files = filter(lambda x: x not in _file_names, input_files.values())
    if len(missing_files):
      for f in missing_files:
        logger.error('File {:s} missing on platform project'.format(f['file_name']))
      return False
    return True





  def job_status(self, job_id):
    """

    :param job_id:
    :return: 'running', 'finished' or 'error'
    """
    if self.job_id[job_id].is_alive():
      return 'running'
    else:
      return 'finished'

  @staticmethod
  def test_process(app, duration, output_files):
    time.sleep(duration)
    for k, v in output_files.items():
      _k = app.get('output_mapping', {}).get(k, None) or k
      open(v, 'w').write('File from {:s} to be named {:s}'.format(_k, v))


