"""
methods for using conda envs from a notebook

simply importing this module updates os.environ['PATH'] to include the current conda env's bin dir.

    from jme.jupy_tools import conda

to switch environments (adjusting sys.path and os.environ['PATH'] use conda.activate

    conda.activate('base')

The argument can be an environment name or path, just as with the command line

to revert, use deactivate:

    conda.deactivate()

Notes:

 * this module remembers all the past activate() calls and deactivate() steps back through them
 * activate(env) does not change the python module path by default, but it can (see docs)

"""

import sys, os, re
import subprocess

## Start with a hack for a special case that I come across often:
#    module being imported to a jupyter env that is run in a custom conda environment
#    BUT: the shell path doesn't reflect the env
#  So:
#   Figure out what env it is from the python path
#   Update the shell path

# these will be set at the end of loading, but should be constants after that
CONDA_BIN_DIR = None
ENVS = {}

# save the shell path that's not part of conda or jupyter
NON_CONDA_SHELL_PATH = [
    p
    for p in os.environ["PATH"].split(":")
    if re.search(r"(jupyter|conda)", p, re.I) is None
]
# save the original python path
ORIG_PATH = sys.path

# keep track of (de)activations
PREVIOUS_ENVS = []
active_env = None

MODULE_PATH = os.path.dirname(os.path.abspath(__file__))


# Functions for inspecting the local conda setup
def get_conda_bin_path():
    # is conda in the path?
    try:
        conda_path = subprocess.check_output(['which', 'conda'], ).rstrip()
    except CalledProcessError:
        # conda not installed in this shell
        return None
    
    try:
        # subprocess usually returns a byte string, decode it
        conda_path = conda_path.decode()
    except AttributeError:
        # It was already a str
        pass
    
    # we want the bin dir, not the full path to the conda executable
    return os.path.dirname(conda_path)

def discover_named_envs():
    env_bytes = subprocess.check_output(['conda', 'info', '--envs'])
    if isinstance(env_bytes, bytes):
        env_bytes = env_bytes.decode()
    for line in env_bytes.split("\n"):
        # skip comments
        if line.startswith("#"):
            continue
        if len(line.strip()) == 0:
            continue
        m = re.search(
            r'^(\S*)\s+\*?\s+(/.+)$',
            line
        )
        if m:
            name, path = m.groups()
            if len(name) > 0:
                ENVS[name] = path.strip()
            else:
                name = guess_env_name(path)
                if name is not None and name not in ENVS:
                    ENVS[name] = path

def guess_env_name(env_path):
    dname, bname = os.path.split(env_path)
    # skip snakemake envs
    if re.search('snakemake', dname, re.I):
        return None
    # if "conda" is in the name of the parent folder
    if re.search('conda', os.path.basename(dname), re.I):
        # use the base name
        return bname
    return None


def list_envs():
    return list(ENVS)

def get_active_env():
    return active_env

def activate(env, set_python_path=False, set_shell_path=True, debug=False):
    """
    This "activates" a requested conda environment. By default, the activation
    is only partial and only affects the shell's PATH environment variable
    (for more see set_python_path and set_shell_path params below).

    This let's you easily access command line programs from other conda environments.

    params:
     * env: the name or path of a conda environment. As of now, names are not perfectly
       supported. Absolute paths should always work.
     * set_python_path: by default, pythons sys.path is NOT modified. To modify sys.path
       to be able to load python modules from this env, pass set_python_path=True.
       (This can cause problems, particularly if you switch envs frequently)
     * set_shell_path: by default, the os.environ[PATH] is updated to include the env's
       bin folder. Set this to False to leave the PATH unchanged.

    """
    global active_env
    current_env = {"env": active_env, "path": os.environ["PATH"], "pypath": sys.path}

    # get the env location
    try:
        # is it a named env?
        env = ENVS[env]
    except KeyError:
        if not os.path.isdir(f"{env}/bin"):
            # if it's a path, it must have a bin dir
            raise Exception(f"Cannot find env: {env}")
    # add notebook's base dir as a fall back for shell PATH
    env_dirs = [env, ]
    if CONDA_BIN_DIR is not None:
        env_dirs.append(CONDA_BIN_DIR)

    # save env path to global var
    active_env = env

    # add /bin to end of env dirs (and use aboslute paths)
    env_dirs = [os.path.abspath(d) + "/bin" for d in env_dirs]

    # update os.environ PATH (prepend env's bin dir)
    if set_shell_path:
        os.environ["PATH"] = ":".join(env_dirs + NON_CONDA_SHELL_PATH)

    # update python's sys.path
    if set_python_path:
        if env is None:
            sys.path = ORIG_PATH
        else:
            # get sys.path from env's Python executable
            cmd = env_dirs[0] + '/python -c "import sys; print(sys.path)"'
            sys.path = eval(subprocess.check_output(cmd, shell=True))
            # add path to this module
            sys.path.append(MODULE_PATH)

    # save previous state
    PREVIOUS_ENVS.append(current_env)


def deactivate():
    """ revert to previous state of sys.path and os.environ[PATH] """
    prev_env = PREVIOUS_ENVS.pop(-1)
    os.environ["PATH"] = prev_env["path"]
    sys.path = prev_env["pypath"]
    active_env = prev_env['env']

def fix_env(debug=False):
    """ attempt to add the current env to the os.environ(PATH) """

    # look for things of the form /some/path/lib/python...
    #  ...where there's a coda-meta dir in the same dir as /lib/
    for path in sys.path:
        if debug:
            print(f'Found "{path}" in sys.path')
        m = re.search(r'(.+)/(lib/python|lib_pypy|site-packages).*', path)
        if m:
            prefix = m.group(1)
            conda_meta = os.path.join(prefix, 'conda-meta')
            if os.path.exists(conda_meta) and os.path.isdir(conda_meta):
                env = prefix
                break
    else:
        raise Exception("Couldn't find anything in sys.path that looks like conda.")
    
    if debug:
        print(f'activating: {env}')
    activate(env, debug=True)

    
# Attempt to do some things automatically, but don't fail!
try:
    # get the location of conda
    CONDA_BIN_DIR = get_conda_bin_path()

    # get the list of installed envs
    discover_named_envs()
    
    # always try to fix the env
    fix_env()
except:
    # but don't fail
    print("WARNING: couldn't auto-acticate current kernel env")
