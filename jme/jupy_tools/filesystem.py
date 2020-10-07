"""
FileSystem related utilities:

 find(path, filters=None):
    recursively returns every file in path that passes given filters
    filters is a list of callables that return True/False based on a path

 glob_wildcards(template, constraints=None):
    finds all paths that match given template
    yields path, wildcards pairs where wildcards is a dict of the values that will furn the template into path using template.format(**wildcards)

    example:
        > for path, wildcards in filesystem.glob_wildcards('tests/{mod}_test.py'):
        -     print(f"mod={wildcards['mod']}: {path}")
        mod=fs: tests/fs_test.py
        mod=cdhit: tests/cdhit_test.py
        mod=utils: tests/utils_test.py

"""

import os, re, glob
from collections import namedtuple

def find(path, filters=None):
    """ recursively find files that match executable filters """
    for root, dirs, files in os.walk(path, topdown=False):
        for name in files:
            if _check_path(root, name, filters):
                yield os.path.join(root, name)

        skips = []
        for name in dirs:
            if _check_path(root, name, filters):
                yield os.path.join(root, name)
            else:
                skips.append(name)

        for name in skips:
            dirs.remove(name)

def get_rexp_filter(rexp):
    if isinstance(rexp, str):
        rexp = re.compile(rexp)
    def check_path_against_rexp(root, name):
        return rexp.search(os.path.join(root, name)) is not None
    return check_path_against_rexp

def _check_path(root, name, filters=None):
    if filters:
        for test in filters:
            if not test(root, name):
                return False
    return True


# matches the bracketed names in a wildcard glob string
TEMPLATE_REXP = re.compile(r"(?<!{)\{\s*(\w+)\s*\}(?!})")


def glob_wildcards(template, constraints=None, as_tuple=False):
    """ should work like the snakemake function:
          * given a template like path/{variable}.ext
          * find the values for variable that match files
          
        Except the return value is different.
        
        Generates for each found file:
            (file_path, dict(wildcard values))
            
        If as_tuple is True, generates named tuple instead of dict.
            
    """
    if constraints is None:
        constraints = {}

    # simple glob for finding files {xxx} -> *
    glob_string = TEMPLATE_REXP.sub("*", template)

    def wc_repl(match):
        """ replace {xxx} with named regex pattern using any constraints """
        wc_name = match.group(1)
        wc_patt = constraints.get(wc_name, r".+")
        return f"(?P<{wc_name}>{wc_patt})"

    # regex for getting wildcards from path
    wildcard_rexp = re.compile(TEMPLATE_REXP.sub(wc_repl, _hide_dots(template)))

    # create named tuple class for returned data
    if as_tuple:
        Wildcards = namedtuple(
            "Wildcards", [m.group(1) for m in TEMPLATE_REXP.finditer(template)]
        )

    # loop over matched files
    for glob_file in glob.glob(glob_string):
        m = wildcard_rexp.match(glob_file)
        if m:
            # transform dict of names->matches to named tuple, if asked
            wildcards = Wildcards(**m.groupdict()) if as_tuple else m.groupdict()
            yield glob_file, wildcards
        else:
            print("WARNING: {} doesn match {}".format(glob_file, wildcard_rexp.pattern))


def _hide_dots(path):
    return re.sub(r"\.", "\.", path)


from numpy import log, power, abs

LN_BASE = log(power(1024, 1 / 3))


def human_readable_bytes(byt):
    """ fixed version of https://stackoverflow.com/a/17754143/663466
     hybrid of https://stackoverflow.com/a/10171475/2595465
      with https://stackoverflow.com/a/5414105/2595465  """
    # return bytes if small
    if byt <= 99:
        return str(int(byt))
    magnitude = int(log(abs(byt)) / LN_BASE)
    if magnitude > 19:
        float_fmt = "%i"
        illion = 20 // 3
    else:
        mag3 = (magnitude + 1) % 3
        float_fmt = "%" + str(mag3) + "." + str(3 - mag3) + "f"
        illion = (magnitude + 1) // 3
    format_str = float_fmt + ["", "K", "M", "G", "T", "P", "E"][illion]
    return (format_str % (byt * 1.0 / (1024 ** illion))).lstrip("0")
