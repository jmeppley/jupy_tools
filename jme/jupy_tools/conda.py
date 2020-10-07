"""
methods for invoking conda from a notebook

simply importing this module updates os.environ['PATH'] to include the current conda env's bin dir.

    from jme.jupy_tools import conda

to switch environments (adjusting sys.path and os.environ['PATH'] use conda.activate

    conda.activate('base')

The argument can be an environment name or path, just as with the command line

to revert, use deactivate:

    conda.deactivate()

"""
# remove anything conda or jupyter related from the path, and add conda from my lysine folder
#  then save to the kernel's env
import sys, os, re
import subprocess

CONDA_BASE_DIR = "/mnt/lysine/jmeppley/miniconda"
OTHER_CONDA_DIRS = [
    "/data/jmeppley/miniconda3",
]
ENVS = {}
for base_dir in [CONDA_BASE_DIR] + OTHER_CONDA_DIRS:
    for d in os.listdir(f"{base_dir}/envs"):
        if os.path.isdir(f"{base_dir}/envs/{d}/bin"):
            ENVS[d] = f"{base_dir}/envs/{d}"
NON_CONDA_PATH = [
    p
    for p in os.environ["PATH"].split(":")
    if re.search(r"(jupyter|conda)", p, re.I) is None
]
ORIG_PATH = sys.path
PREVIOUS_ENVS = []
MODULE_PATH = os.path.dirname(os.path.abspath(__file__))


def list_envs():
    return list(ENVS)


def activate(env=None):
    current_env = {"path": os.environ["PATH"], "pypath": sys.path}

    # update PATH
    if env != None:
        try:
            # is it a named env?
            env = ENVS[env]
        except KeyError:
            if not os.path.isdir(f"{env}/bin"):
                # if it's a path, it must have a bin dir
                raise Exception(f"Cannot find env: {env}")
        env_dirs = [env, CONDA_BASE_DIR]
    else:
        env_dirs = [CONDA_BASE_DIR]

    # add /bin to end of env dirs
    env_dirs = [d + "/bin" for d in env_dirs]

    # set new PATH
    os.environ["PATH"] = ":".join(env_dirs + NON_CONDA_PATH)

    # update sys.path
    if env is None:
        sys.path = ORIG_PATH
    else:
        # get sys.path from Python executable
        cmd = 'python -c "import sys; print(sys.path)"'
        sys.path = eval(subprocess.check_output(cmd, shell=True))
        # add path to this module
        sys.path.append(MODULE_PATH)

    # save previous state
    PREVIOUS_ENVS.append(current_env)


def deactivate():
    prev_env = PREVIOUS_ENVS.pop(-1)
    os.environ["PATH"] = prev_env["path"]
    sys.path = prev_env["pypath"]

def fix_env():
    """ attempt to add the current env to the os.environ(PATH) """
    
    # look for things of the form /some/path/lib/python...
    #  ...where there's a coda-meta dir in the same dir as /lib/
    for path in sys.path:
        halves = path.split('/lib/python', 1)
        if len(halves) > 1:
            half = halves[0]
            conda_meta = os.path.join(half, 'conda-meta')
            if os.path.exists(conda_meta) and os.path.isdir(conda_meta):
                env = half
                break
    else:
        raise Exception("Couldn't find anything in sys.path that looks like conda.")
    activate(half)
    
# always try to fix the env
try:
    fix_env()
except:
    # but don't fail
    print("WARNING: couldn't auto-acticate current kernel env")
