"""
Class that will inspect a jupyter server (via the API) and the local running processes
(via subprocess and `top`) trying to connect kernel processes to notebook paths
"""

import datetime
import os
import re
import subprocess

import pandas
import requests

MEM_KEY = '%MEM'
RES_KEY = 'RES'
CPU_KEY = '%CPU'
TOP_KEYS = [CPU_KEY, RES_KEY, MEM_KEY]

TOP_CMD_TEMPLATE = 'top -b -n 1 -u {user} -o {sort_key} -c -w 333 | grep -m {N} json$'
TOP_HEADER_CMD = 'top -b -n 1 -c -w 333 | grep -m 1 COMMAND$'

def decode(possible_bytes):
    """ decode bytes to string. Don't complain if it's already a string """
    try:
        return possible_bytes.decode()
    except AttributeError:
        # subprocess returned a string, bot bbytes. OK
        return possible_bytes

# inspect the top header line to figure out where the MEM and CPU values are in a line
def get_cpu_mem_posns_in_top():
    top_header = decode(subprocess.check_output(TOP_HEADER_CMD, shell=True)).rstrip()

    # this only works for right justified (numeric?) fields
    return {key:re.search(f"\s+{key}", top_header).span() for key in TOP_KEYS}
  

def get_top_kernels_from_top(user, sort_key, N):
    """ returns the first N top entries sorted by sort_key """
    output = subprocess.check_output(
        TOP_CMD_TEMPLATE.format(user=user, sort_key=sort_key, N=N),
        shell=True
    )
    return decode(output).rstrip().split('\n')


class KernelFinder():
    def __init__(self, 
                 jupyter_url,
                 jupyter_token,
                 notebook_root=os.environ['HOME'],
                 username=os.getlogin()):
        self.url = jupyter_url
        self.token = jupyter_token
        self.user = username
        self.root = notebook_root
        self.sessions_by_kernel = self.get_sessions_by_kernel(use_cache=False)
        self.col_spans = self.get_col_spans(use_cache=False)

    def get_kernels_table(self, **kwargs):
        """
        returns a DataFrame of kernels (1-per-row) with columns:
         * hash
         * notebook path
         * notebook mtime
         * RAM usage (%)
         * RAM usage (#)
         * CPU usage (%)
        
        Output can be filtered or sorted as in det_top_N_kernels
        """
        some_list = [1,2,3,4]
        data_generator = (
            (kernel, 
             kdata['path'], 
             self.parse_vals_from_top_line(kdata['top_line'],
                                           use_cache=kwargs.get('use_cacke', True)),
            )
            for kernel, kdata
            in self.get_top_N_kernels(**kwargs).items()
        )
        row_generator = (
            [kernel, *os.path.split(nb_path), self._getmtime(nb_path),] + [top_data[k] for k in TOP_KEYS]
            for kernel, nb_path, top_data
            in data_generator
        )
        df = pandas.DataFrame(
            row_generator,
            columns=(['kernel', 'nb_path', 'nb_name', 'nb_mtime'] + TOP_KEYS),
        ).set_index('kernel')
        df['age'] = datetime.datetime.now() - df.nb_mtime
        df[CPU_KEY] = [float(c) for c in df[CPU_KEY]]
        df[MEM_KEY] = [float(m) for m in df[MEM_KEY]]
        return df
    
    def _getmtime(self, notebook_path):
        try:
            return pandas.to_datetime(
                datetime.datetime.fromtimestamp(
                    os.path.getmtime(os.path.join(self.root, notebook_path))
                )
            )
        except FileNotFoundError:
            return None
        
    def get_top_N_kernels(self, N=-1, key=MEM_KEY, use_cache=True):
        """
        Get the top N kernels by CPU or MEM (defaults to all (N = -1))
        
        if use_cache is False, re-query the jupyter server for the latest list of kernels
        """        
        
        session_by_kernel = self.get_sessions_by_kernel(use_cache)

        # !top -b -n 1 -u jmeppley -o {key} -c -w 333 | grep -m {N} json$
        top_N = get_top_kernels_from_top(self.user, key, N)

        kernels = {}
        for line in top_N:
            try:
                kernel = re.search(r'runtime/kernel-(\S+).json', line).group(1)
            except AttributeError:
                print(f"WARNING: can't parse line: {line}")
                continue
            try:
                kernel_path = session_by_kernel[kernel]['notebook']['path']
            except KeyError:
                kernel_path = "**missing kernel**"
            kernels[kernel] = {'path': kernel_path, 'top_line': line}
        return kernels

    def get_sessions_by_kernel(self, use_cache=True):
        if use_cache is False:
            response = requests.get(f"{self.url}/api/sessions", 
                                    verify=False,
                                    headers={"Authorization": f"token {self.token}"})
            self.sessions_by_kernel = {
                s['kernel']['id']: s
                for s in response.json()
            }
        return self.sessions_by_kernel
    
    def get_col_spans(self, use_cache=True):
        if use_cache is False:
            self.col_spans = get_cpu_mem_posns_in_top()
        return self.col_spans

    def parse_vals_from_top_line(self, top_line, use_cache=True):
        col_spans = self.get_col_spans(use_cache)
        
        return {
            key:top_line[span[0]:span[1]].strip()
            for key, span in col_spans.items()
        } 

