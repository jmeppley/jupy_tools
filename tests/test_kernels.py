"""
kernels.py talks to a jupyter server (with requests) and runs `top`. These
tests replace both with canned output.
"""
import datetime
import importlib
import os
import subprocess
import sys

import pandas
import pytest


@pytest.fixture
def kernels():
    from jme.jupy_tools import kernels
    return kernels


def test_import_without_login(monkeypatch):
    """ importing doesn't need a controlling terminal (eg in CI) """
    def no_login():
        raise OSError(25, "Inappropriate ioctl for device")
    monkeypatch.setattr(os, 'getlogin', no_login)
    monkeypatch.delitem(sys.modules, 'jme.jupy_tools.kernels', raising=False)
    try:
        importlib.import_module('jme.jupy_tools.kernels')
    finally:
        sys.modules.pop('jme.jupy_tools.kernels', None)


def test_import_without_home(monkeypatch):
    monkeypatch.delenv('HOME', raising=False)
    monkeypatch.delitem(sys.modules, 'jme.jupy_tools.kernels', raising=False)
    try:
        importlib.import_module('jme.jupy_tools.kernels')
    finally:
        sys.modules.pop('jme.jupy_tools.kernels', None)


def test_decode(kernels):
    assert kernels.decode(b"abc\n") == "abc\n"
    assert kernels.decode("abc") == "abc"


# canned `top -b -c -w 333` output (right justified numeric columns)
TOP_HEADER = "    PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND"
TOP_LINE = "{pid:>7} tester    20   0 {virt:>7} {res:>6} {shr:>6} S {cpu:>5} {mem:>5}   1:23.45 /opt/conda/bin/python -m ipykernel_launcher -f /home/tester/.local/share/jupyter/runtime/kernel-{kernel}.json"
TOP_LINES = [
    TOP_LINE.format(pid=101, virt='5.0g', res='2.1g', shr='10m', cpu='99.0',
                    mem='12.5', kernel='aaaa-1111'),
    TOP_LINE.format(pid=102, virt='1.0g', res='300m', shr='10m', cpu='0.0',
                    mem='1.5', kernel='bbbb-2222'),
    TOP_LINE.format(pid=103, virt='1.0g', res='100m', shr='10m', cpu='5.0',
                    mem='0.5', kernel='cccc-3333'),
]


@pytest.fixture
def fake_top(kernels, monkeypatch):
    calls = []

    def check_output(cmd, shell=False):
        calls.append(cmd)
        if cmd == kernels.TOP_HEADER_CMD:
            return (TOP_HEADER + "\n").encode()
        return ("\n".join(TOP_LINES) + "\n").encode()

    monkeypatch.setattr(subprocess, 'check_output', check_output)
    return calls


def test_get_cpu_mem_posns_in_top(kernels, fake_top):
    spans = kernels.get_cpu_mem_posns_in_top()
    assert set(spans) == set(kernels.TOP_KEYS)
    for key, (start, end) in spans.items():
        assert TOP_HEADER[start:end].strip() == key
        # the same span pulls the value out of a data line
        value = TOP_LINES[0][start:end].strip()
        assert value == {'%CPU': '99.0', '%MEM': '12.5', 'RES': '2.1g'}[key]


def test_get_top_kernels_from_top(kernels, fake_top):
    lines = kernels.get_top_kernels_from_top('tester', kernels.MEM_KEY, 3)
    assert lines == TOP_LINES
    assert fake_top == [kernels.TOP_CMD_TEMPLATE.format(
        user='tester', sort_key='%MEM', N=3)]
    assert fake_top[0] == \
        'top -b -n 1 -u tester -o %MEM -c -w 333 | grep -m 3 json$'


class FakeResponse:
    def __init__(self, data):
        self.data = data

    def json(self):
        return self.data


SESSIONS = [
    {'kernel': {'id': 'aaaa-1111'}, 'notebook': {'path': 'proj/big.ipynb'}},
    {'kernel': {'id': 'bbbb-2222'}, 'notebook': {'path': 'small.ipynb'}},
    {'kernel': {'id': 'dddd-4444'}, 'notebook': {'path': 'idle.ipynb'}},
]


@pytest.fixture
def finder(kernels, fake_top, monkeypatch, tmp_path):
    requests_made = []

    def fake_get(url, verify=True, headers=None):
        requests_made.append((url, headers))
        return FakeResponse(SESSIONS)

    monkeypatch.setattr(kernels.requests, 'get', fake_get)
    (tmp_path / 'proj').mkdir()
    (tmp_path / 'proj' / 'big.ipynb').write_text('{}')
    (tmp_path / 'small.ipynb').write_text('{}')
    finder = kernels.KernelFinder('https://jupyter:8888', 'TOKEN',
                                  notebook_root=str(tmp_path),
                                  username='tester')
    finder.requests_made = requests_made
    return finder


def test_kernel_finder_defaults(kernels, fake_top, monkeypatch, tmp_path):
    """ the user and notebook root default to the current user and home dir """
    monkeypatch.setattr(kernels.requests, 'get',
                        lambda url, verify=True, headers=None:
                        FakeResponse(SESSIONS))

    def no_login():
        raise OSError(25, "Inappropriate ioctl for device")
    monkeypatch.setattr(os, 'getlogin', no_login)
    # getpass.getuser() checks these in order
    for var in ['LOGNAME', 'USER', 'LNAME', 'USERNAME']:
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv('LOGNAME', 'env_user')

    # looked up when the finder is made, not when the module was imported
    monkeypatch.setenv('HOME', str(tmp_path / 'home1'))
    finder = kernels.KernelFinder('https://jupyter:8888', 'TOKEN')
    assert finder.user == 'env_user'
    assert finder.root == str(tmp_path / 'home1')

    monkeypatch.setenv('HOME', str(tmp_path / 'home2'))
    monkeypatch.setenv('LOGNAME', 'other_user')
    finder = kernels.KernelFinder('https://jupyter:8888', 'TOKEN')
    assert finder.user == 'other_user'
    assert finder.root == str(tmp_path / 'home2')


def test_kernel_finder_init(finder):
    assert finder.user == 'tester'
    assert finder.requests_made == [
        ('https://jupyter:8888/api/sessions', {'Authorization': 'token TOKEN'})]
    assert set(finder.sessions_by_kernel) == {'aaaa-1111', 'bbbb-2222',
                                              'dddd-4444'}
    assert set(finder.col_spans) == {'%CPU', '%MEM', 'RES'}


def test_kernel_finder_caching(finder):
    finder.get_sessions_by_kernel()
    assert len(finder.requests_made) == 1
    finder.get_sessions_by_kernel(use_cache=False)
    assert len(finder.requests_made) == 2


def test_parse_vals_from_top_line(finder):
    assert finder.parse_vals_from_top_line(TOP_LINES[1]) == \
        {'%CPU': '0.0', '%MEM': '1.5', 'RES': '300m'}


def test_get_top_n_kernels(finder):
    found = finder.get_top_N_kernels()
    assert list(found) == ['aaaa-1111', 'bbbb-2222', 'cccc-3333']
    assert found['aaaa-1111']['path'] == 'proj/big.ipynb'
    assert found['aaaa-1111']['top_line'] == TOP_LINES[0]
    # running kernel without a session
    assert found['cccc-3333']['path'] == '**missing kernel**'


def test_get_top_n_kernels_bad_line(finder, monkeypatch, capsys, kernels):
    monkeypatch.setattr(kernels, 'get_top_kernels_from_top',
                        lambda user, key, N: ['garbage line', TOP_LINES[0]])
    found = finder.get_top_N_kernels()
    assert list(found) == ['aaaa-1111']
    assert "can't parse line: garbage line" in capsys.readouterr().out


def test_getmtime(finder, tmp_path):
    mtime = finder._getmtime('small.ipynb')
    expected = datetime.datetime.fromtimestamp(
        os.path.getmtime(tmp_path / 'small.ipynb'))
    assert mtime == pandas.Timestamp(expected)
    assert finder._getmtime('no_such.ipynb') is None


def test_get_kernels_table(finder):
    table = finder.get_kernels_table()
    assert list(table.index) == ['aaaa-1111', 'bbbb-2222', 'cccc-3333']
    assert list(table.columns) == ['nb_path', 'nb_name', 'nb_mtime', '%CPU',
                                   'RES', '%MEM', 'age']
    assert table.loc['aaaa-1111', 'nb_path'] == 'proj'
    assert table.loc['aaaa-1111', 'nb_name'] == 'big.ipynb'
    assert table.loc['bbbb-2222', 'nb_path'] == ''
    assert table.loc['aaaa-1111', '%CPU'] == 99.0
    assert table.loc['bbbb-2222', '%MEM'] == 1.5
    assert table.loc['aaaa-1111', 'RES'] == '2.1g'
    assert table.loc['aaaa-1111', 'age'] >= pandas.Timedelta(0)
    # no notebook file for the missing kernel
    assert pandas.isna(table.loc['cccc-3333', 'nb_mtime'])
