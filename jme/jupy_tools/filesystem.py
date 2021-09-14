"""
FileSystem related utilities:

  find(path, filters=None):
     recursively returns every file in path that passes given filters
     filters is a list of callables that return True/False based on a path

  glob_wildcards(template, constraints=None):
     finds all paths that match given template
     yields (path, wildcards) pairs
        where wildcards is a dict of the values that will furn the template into path
        using template.format(**wildcards)

     example:
        > for path, wildcards in filesystem.glob_wildcards('tests/{mod}_test.py'):
        -     print(f"mod={wildcards['mod']}: {path}")
        mod=fs: tests/fs_test.py
        mod=cdhit: tests/cdhit_test.py
        mod=utils: tests/utils_test.py

  human_readable_bytes(byt):
      return a string like 5.6G or 1.0K or 45T given a number of bytes

  get_file_sizes_and_dates_by_uid(volume, users=None, min_age=0, follow_links=False):
      walk a directory tree to catalog contents by size/date/user

  get_file_size_table(usage_data, min_date):
      aggregate total file sizes by date and user
"""

import os, re, glob, pandas, numpy
from datetime import datetime
from collections import namedtuple, defaultdict

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


def glob_wildcards(template, constraints=None, as_tuple=False, debug=False):
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

    named_patterns = set()
    def wc_repl(match):
        """ replace {xxx} with named regex pattern using any constraints """
        wc_name = match.group(1)
        if wc_name in named_patterns:
            return f"(?P={wc_name})"
        named_patterns.add(wc_name)
        wc_patt = constraints.get(wc_name, r".+")
        return f"(?P<{wc_name}>{wc_patt})"

    # regex for getting wildcards from path
    wildcard_pattern = TEMPLATE_REXP.sub(wc_repl, _hide_dots(template))
    if debug:
        print(f"Wildcard regexp: '{wildcard_pattern}'")
    wildcard_rexp = re.compile(wildcard_pattern)
    
    # create named tuple class for returned data
    if as_tuple:
        Wildcards = namedtuple(
            "Wildcards", list(set(m.group(1) for m in TEMPLATE_REXP.finditer(template)))
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


def get_file_sizes_and_dates_by_uid(volume, users=None, min_age=0, follow_links=False, time_stat='mtime'):
    """ Collect date and size by user id 
    
    Params:
    
     * users: only consider these users (ids or names)
     * min_age: only save files at least this old (in seconds)
     * follow_links: passed to os.walk (default is False)
     * time_stat: one of 'mtime', 'ctime', 'atime', or 'max'
    """

    # translate user ids to names
    userid_map = get_user_lookup_table().to_dict()

    # translate userids to names in include list
    if users is not None:
        users = set(userid_map.get(u, u) for u in users)

    usage_data = defaultdict(lambda: [])
    min_date = int(datetime.now().timestamp())
    now = datetime.now().timestamp()
    for root_path, folder_list, file_list in os.walk(volume, followlinks=follow_links):
        for file_name in file_list:
            try:
                file_path = os.path.join(root_path, file_name)
                if not(os.path.isfile(file_path)):
                    # skip broken links
                    continue
                file_stats = os.stat(file_path)

                # filter by owner if user list given
                ownerid = file_stats.st_uid
                owner = userid_map.get(ownerid, ownerid)
                if users is not None and owner not in users:
                    continue

                # get the user selected time stat
                mtime = file_stats.st_ctime if time_stat == 'ctime' \
                    else file_stats.st_mtime if time_stat == 'mtime' \
                    else file_stats.st_atime if time_stat == 'atime' \
                    else max(file_stats.st_mtime, file_stats.st_ctime, file_stats.st_atime)
                
                # keep track of oldest file
                min_date = min(mtime, min_date)
                # filter by age
                file_age = now - mtime
                if file_age < min_age:
                    continue
                usage_data[owner].append((file_stats.st_size,
                                          mtime,
                                          file_path,
                                         ))
            except:
                pass

    return usage_data, min_date


TIME_SPANS = {
    'minutes': 60,
    'hours': 3600,
    'days': 3600*24,
    'weeks': 3600*24*7,
    'months': 3600*24*30,
    'years': 3600*24*365,
}

def get_file_size_table(usage_data, min_date,
                        age_bin_size=2,
                        age_bin_type='weeks', min_age=0):
    """ translate files sizes and dates into table """

    now = datetime.now().timestamp()
    if age_bin_type not in TIME_SPANS:
        raise Exception("I don't know the time span {}. Please specify one of: {}".format(
            age_bin_type,
            ", ".join(TIME_SPANS.keys()),
        ))

    age_bins_step = age_bin_size * TIME_SPANS[age_bin_type]
    oldest_age = now - min_date
    age_bin_bounds = numpy.arange(0, oldest_age + age_bins_step, age_bins_step)
    
    counts = {}
    now = datetime.now().timestamp()
    for owner, file_data_list in usage_data.items():
        owner_counts = counts.setdefault(owner, {})
        for file_data in file_data_list:
            size = file_data[0]
            file_age = now - file_data[1]
            if file_age < min_age:
                continue
            age_bin = int(file_age/age_bins_step)
            owner_counts[age_bin] = owner_counts.get(age_bin, 0) + size
            
    
    # make into a data frame
    file_size_table = pandas.DataFrame(counts)
    # headers...
    #users = get_user_lookup_table()
    #file_size_table.columns = [users.get(c,c) for c in file_size_table.columns]
    
    file_size_table.index = \
            [get_bin_bounds_string(i, 
                                   age_bin_bounds, 
                                   lambda b: \
                                       str(int(b/TIME_SPANS[age_bin_type])), 
                                   "{} old".format(age_bin_type)) \
             for i in file_size_table.index]
    return file_size_table    
    

def get_bin_bounds_string(bin_index, bin_bounds, to_str=repr, suffix=""):
    """ retuns, for example: '15 to 20 months' given:
      bin_index: the location in bin_bounds to find the start value
      bin_bounds: a list of bounding values
      suffix: the class of the bounds. EG 'months' or 'days'
    """
    return "{} to {} {}".format(to_str(bin_bounds[bin_index]), to_str(bin_bounds[bin_index + 1]), suffix)

def get_user_lookup_table():
    """ returns series mapping user id to user name """
    return pandas.read_table('/etc/passwd', sep=':', names=['user','code','id','group','home','shell'], index_col=2)['user']
