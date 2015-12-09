"""
To use this set the CGC_TOKEN environment variable
"""

import os
import json
import time
from collections import OrderedDict
import logging
import warnings

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
  logger.debug('Path: {:s}\nQuery: {:s}\nData: {:s}\nResponse: {:s}'.format(path, str(query), str(data), str(response_dict)))

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


def create_new_task(project, app, app_list, input_files, output_files, file_list, task_dict, run=True, auth_token=None):
  """

  :param project:
  :param app:
  :param app_list:
  :param input_files:
  :param output_files:
  :param file_list:
  :param task_dict:
  :param run:
  :param auth_token:
  :return:
  """
  def _resolve_input_pin(_k, _app):
    # If it exists, use the input mapping to transform k, otherwise, leave unchanged
    return _app.get('input_mapping', {_k: _k})[_k]

  def _get_file_dict(_input_file, _file_list):
    plat_file = filter(lambda x: x['name'] == _input_file, _file_list)
    if len(plat_file) == 0:
      logger.error('File {:s} missing on platform project'.format(_input_file))
      raise RuntimeError('File {:s} missing on platform project'.format(_input_file))
    plat_file = plat_file[0]
    return {"class": "File", "path": plat_file['id'], "name": plat_file["name"]}

  plat_app = filter(lambda x: x['name'] == app['app_name'], app_list)[0]
  # This should not fail as we have checked that the app exists

  data = {
    "description": "A benchmarking task",
    "name": bench.create_filename_prefix_from_metadata(task_dict['metadata'], use_hash=False),
    "app_id": plat_app['id'],
    "project": project['id'],
    "inputs": {_resolve_input_pin(k, app): _get_file_dict(i_file, file_list) for k, i_file in input_files.items()}
  }

  return api('/tasks', data=data, query={'action': 'run'} if run else None, auth_token=auth_token)


def check_task_status(task_id, auth_token=None):
  status_mapping = {
    'Completed': 'finished',
    'Active': 'running',
    'Aborted': 'error'
  }
  status = api('/tasks/{:s}'.format(task_id))
  return status_mapping[status['status']]


class SBGSDK2Executor(bench.BaseExecutor):
  """Execute benchmarking tasks on an SDK2 enabled platform"""

  def __init__(self, bench_run, project_name, auth_token=None):
    """Setup communication with the platform.

    :param bench_run:
    :param project_name:
    :param auth_token: Platform authentication token
    :return:
    """
    self.auth_token, self.bench_run, self.project, self.project_name = \
      auth_token, bench_run, get_existing_project(project_name, auth_token), project_name
    self._check_project()

    self.apps = get_all_apps_in_project(self.project, self.auth_token)
    self.files_in_project = get_all_files_in_project(self.project, self.auth_token)

    # Check to see all resources are available
    run_name_ok, files_ok, apps_ok = self._check_run_name(), self._check_initial_files(), self._check_apps()
    if not (run_name_ok and files_ok and apps_ok):
      warnings.warn('There were some problems setting up the benchmark. Please see log')

  def _check_project(self):
    """Make sure project exists"""
    if self.project is None:
      logger.error('No project named {:s} found on platform'.format(self.project_name))
      raise RuntimeError('No project named {:s} found on platform'.format(self.project_name))

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

  def start_job(self, app, input_files, output_files, task_dict):
    """Find the app on the platform and start a new task with the given input files.

    :param app:
    :param input_files:
    :param output_files:
    :param task_dict:
    :return:
    """
    # Refresh the file list
    self.files_in_project = get_all_files_in_project(self.project, self.auth_token)
    p = create_new_task(self.project, app, self.apps,
                        input_files, output_files, self.files_in_project,
                        task_dict, auth_token=self.auth_token)
    pid = p['id']
    self.job_id[pid] = p
    self.job_files[pid] = {
      'app': app,
      'input_files': input_files,
      'output_files': output_files
    }
    return pid

  def job_status(self, job_id):
    """

    :param job_id:
    :return: 'running', 'finished' or 'error'
    """
    return check_task_status(job_id, auth_token=self.auth_token)