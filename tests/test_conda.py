import os
import re
import subprocess
import sys

import pytest

from jme.jupy_tools import conda


@pytest.fixture
def clean_conda(monkeypatch):
    """
    Protect sys.path, os.environ['PATH'] and the conda module's globals from
    changes made by the tests.
    """
    monkeypatch.setattr(sys, 'path', list(sys.path))
    monkeypatch.setenv('PATH', os.environ['PATH'])
    monkeypatch.setattr(conda, 'active_env', None)
    monkeypatch.setattr(conda, 'PREVIOUS_ENVS', [])
    monkeypatch.setattr(conda, 'ENVS', {})
    monkeypatch.setattr(conda, 'CONDA_BIN_DIR', None)
    monkeypatch.setattr(conda, 'NON_CONDA_SHELL_PATH', ['/usr/bin', '/bin'])
    return conda


def make_fake_env(root, name, python_path=None):
    """ make a directory that looks enough like a conda env for activate() """
    env = root / name
    (env / 'bin').mkdir(parents=True)
    (env / 'conda-meta').mkdir()
    site_packages = env / 'lib' / 'python3.99' / 'site-packages'
    site_packages.mkdir(parents=True)
    if python_path is not None:
        python = env / 'bin' / 'python'
        python.write_text(f"#!/bin/sh\necho \"{python_path!r}\"\n")
        python.chmod(0o755)
    return env


def test_load_conda():
    """ importing never fails, even without conda """
    from jme.jupy_tools import conda  # noqa: F401


def test_non_conda_shell_path():
    for path in conda.NON_CONDA_SHELL_PATH:
        assert re.search(r'(jupyter|conda)', path, re.I) is None


def test_activate_env():
    """ activate and deactivate a real conda env (if there are any) """
    if len(conda.ENVS) == 0:
        pytest.skip("no conda environments found")
    env_names = list(conda.ENVS)
    env_name = env_names[int(len(conda.ENVS) / 2)]

    base_env = conda.active_env
    orig_path = os.environ['PATH']
    try:
        conda.activate(env_name)
        assert conda.ENVS[env_name] == conda.active_env
    finally:
        conda.deactivate()

    assert base_env == conda.active_env
    assert os.environ['PATH'] == orig_path


###
# inspecting the conda installation


def test_get_conda_bin_path(monkeypatch):
    monkeypatch.setattr(subprocess, 'check_output',
                        lambda cmd: b'/opt/conda/bin/conda\n')
    assert conda.get_conda_bin_path() == '/opt/conda/bin'


def test_get_conda_bin_path_str(monkeypatch):
    monkeypatch.setattr(subprocess, 'check_output',
                        lambda cmd: '/opt/conda/bin/conda')
    assert conda.get_conda_bin_path() == '/opt/conda/bin'


def test_get_conda_bin_path_missing(monkeypatch):
    def no_conda(cmd):
        raise subprocess.CalledProcessError(1, cmd)
    monkeypatch.setattr(subprocess, 'check_output', no_conda)
    assert conda.get_conda_bin_path() is None


CONDA_INFO_ENVS = b"""# conda environments:
#
base                     /opt/conda
analysis              *  /opt/conda/envs/analysis
tools                    /home/me/.conda/envs/tools
                         /home/me/conda_envs/unnamed
                         /home/me/project/.snakemake/conda/abc123
                         /home/me/project/some_env

"""


def test_discover_named_envs(clean_conda, monkeypatch):
    monkeypatch.setattr(subprocess, 'check_output', lambda cmd: CONDA_INFO_ENVS)
    conda.discover_named_envs()
    assert conda.ENVS == {
        'base': '/opt/conda',
        'analysis': '/opt/conda/envs/analysis',
        'tools': '/home/me/.conda/envs/tools',
        'unnamed': '/home/me/conda_envs/unnamed',
    }
    assert sorted(conda.list_envs()) == ['analysis', 'base', 'tools', 'unnamed']


def test_discover_named_envs_no_conda(clean_conda, monkeypatch):
    def no_conda(cmd):
        raise FileNotFoundError(cmd[0])
    monkeypatch.setattr(subprocess, 'check_output', no_conda)
    assert conda.discover_named_envs() is None
    assert conda.ENVS == {}


@pytest.mark.parametrize("path, name", [
    ('/home/me/conda_envs/myenv', 'myenv'),
    ('/home/me/.conda/myenv', 'myenv'),
    ('/home/me/project/.snakemake/conda/abc123', None),
    ('/home/me/project/env', None),
    # the parent folder (not a grandparent) has to look like conda
    ('/opt/miniconda3/envs/myenv', None),
])
def test_guess_env_name(path, name):
    assert conda.guess_env_name(path) == name


###
# activate / deactivate


def test_activate_path(clean_conda, tmp_path):
    env = make_fake_env(tmp_path, 'env1')
    orig_path = os.environ['PATH']
    orig_sys_path = list(sys.path)

    conda.activate(str(env))
    assert conda.get_active_env() == str(env)
    assert os.environ['PATH'] == f"{env}/bin:/usr/bin:/bin"
    # python path is untouched by default
    assert sys.path == orig_sys_path

    conda.deactivate()
    assert conda.get_active_env() is None
    assert os.environ['PATH'] == orig_path
    assert conda.PREVIOUS_ENVS == []


def test_activate_named(clean_conda, tmp_path):
    env = make_fake_env(tmp_path, 'env1')
    conda.ENVS['my_env'] = str(env)
    conda.activate('my_env')
    assert conda.active_env == str(env)
    assert os.environ['PATH'].startswith(f"{env}/bin:")


def test_activate_relative_path(clean_conda, tmp_path, monkeypatch):
    make_fake_env(tmp_path, 'env1')
    monkeypatch.chdir(tmp_path)
    conda.activate('env1')
    assert os.environ['PATH'].startswith(f"{tmp_path}/env1/bin:")


def test_activate_stack(clean_conda, tmp_path):
    """ deactivate() steps back through previous activations """
    env1 = make_fake_env(tmp_path, 'env1')
    env2 = make_fake_env(tmp_path, 'env2')
    orig_path = os.environ['PATH']

    conda.activate(str(env1))
    path1 = os.environ['PATH']
    conda.activate(str(env2))
    assert conda.active_env == str(env2)
    assert len(conda.PREVIOUS_ENVS) == 2

    conda.deactivate()
    assert conda.active_env == str(env1)
    assert os.environ['PATH'] == path1

    conda.deactivate()
    assert conda.active_env is None
    assert os.environ['PATH'] == orig_path

    with pytest.raises(IndexError):
        conda.deactivate()


def test_activate_no_shell_path(clean_conda, tmp_path):
    env = make_fake_env(tmp_path, 'env1')
    orig_path = os.environ['PATH']
    conda.activate(str(env), set_shell_path=False)
    assert conda.active_env == str(env)
    assert os.environ['PATH'] == orig_path


def test_activate_python_path(clean_conda, tmp_path):
    fake_sys_path = ['/fake/lib/python3.99', '/fake/site-packages']
    env = make_fake_env(tmp_path, 'env1', python_path=fake_sys_path)
    orig_sys_path = sys.path

    conda.activate(str(env), set_python_path=True)
    assert sys.path == fake_sys_path + [conda.MODULE_PATH]

    conda.deactivate()
    assert sys.path is orig_sys_path


def test_activate_missing_env(clean_conda, tmp_path):
    with pytest.raises(Exception, match="Cannot find env"):
        conda.activate('no_such_env_name')
    # a dir without a bin/ folder is not an env
    with pytest.raises(Exception, match="Cannot find env"):
        conda.activate(str(tmp_path))
    assert conda.active_env is None
    assert conda.PREVIOUS_ENVS == []


def test_activate_adds_conda_bin_dir(clean_conda, tmp_path):
    """ conda's own bin dir is a fall back after the env's bin dir """
    env = make_fake_env(tmp_path, 'env1')
    conda.CONDA_BIN_DIR = '/opt/conda/bin'
    conda.activate(str(env))
    assert os.environ['PATH'] == f"{env}/bin:/opt/conda/bin:/usr/bin:/bin"

    # deactivating removes it again
    conda.deactivate()
    assert '/opt/conda/bin' not in os.environ['PATH'].split(":")


def test_activate_finds_conda_bin_programs(clean_conda, tmp_path, monkeypatch):
    """
    Programs in the env come first, then programs next to the conda
    executable (CONDA_BIN_DIR comes from `which conda`)
    """
    import shutil
    base_bin = tmp_path / 'miniconda3' / 'bin'
    base_bin.mkdir(parents=True)
    env = make_fake_env(tmp_path, 'env1')
    for path in [base_bin / 'conda', base_bin / 'only_in_base',
                 base_bin / 'in_both', env / 'bin' / 'in_both']:
        path.write_text("#!/bin/sh\n")
        path.chmod(0o755)

    monkeypatch.setattr(subprocess, 'check_output',
                        lambda cmd: f"{base_bin}/conda\n".encode())
    conda.CONDA_BIN_DIR = conda.get_conda_bin_path()
    assert conda.CONDA_BIN_DIR == str(base_bin)

    conda.activate(str(env))
    assert shutil.which('only_in_base') == str(base_bin / 'only_in_base')
    assert shutil.which('in_both') == str(env / 'bin' / 'in_both')
    assert shutil.which('conda') == str(base_bin / 'conda')


def test_activate_relative_conda_bin_dir(clean_conda, tmp_path, monkeypatch):
    env = make_fake_env(tmp_path, 'env1')
    monkeypatch.chdir(tmp_path)
    conda.CONDA_BIN_DIR = 'conda/bin'
    conda.activate(str(env))
    assert os.environ['PATH'].split(":")[1] == f"{tmp_path}/conda/bin"


###
# fix_env


def test_fix_env(clean_conda, tmp_path, capsys):
    env = make_fake_env(tmp_path, 'kernel_env')
    # what sys.path looks like in a kernel running from a conda env
    sys.path[:] = ['/usr/lib/python3.99',
                   f"{env}/lib/python399.zip",
                   f"{env}/lib/python3.99",
                   f"{env}/lib/python3.99/site-packages",
                   '/somewhere/else']
    conda.fix_env(debug=True)
    assert conda.active_env == str(env)
    assert os.environ['PATH'].startswith(f"{env}/bin:")
    out = capsys.readouterr().out
    assert f"activating: {env}" in out


def test_fix_env_no_conda(clean_conda, tmp_path):
    # looks like python, but there is no conda-meta folder
    sys.path[:] = [f"{tmp_path}/lib/python3.99/site-packages"]
    with pytest.raises(Exception, match="Couldn't find anything"):
        conda.fix_env()
