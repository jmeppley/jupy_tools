import sys, os
PY_PATH = sys.path
OS_PATH = os.environ['PATH']

def test_load_conda():
    from jme.jupy_tools import conda

def test_activate_env():
    from jme.jupy_tools import conda

    env_names = list(conda.ENVS)
    env_name = env_names[int(len(conda.ENVS)/2)]

    base_env = conda.active_env
    conda.activate(env_name)
   
    assert conda.ENVS[env_name] == conda.active_env

    conda.deactivate()

    assert base_env == conda.active_env





