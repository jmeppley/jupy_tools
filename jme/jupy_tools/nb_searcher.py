import attr
import os
import re

from . import filesystem

### Constants

HOME = os.environ['HOME']

SEARCH_FILE='file'
SEARCH_CONTENTS='text'

# filters to get just notebooks
_NB_REXP = re.compile(r'\.ipynb$')
_CP_REXP = re.compile(r'\.ipynb_checkpoint')

def _process_filter(filter_param, func_type=SEARCH_FILE):
    """
    standardize the form of the filter parameter.
    
    Take:
     * string (a string or regular expression to search for)
     * re.Pattern (compiled regular expression)
     * a filter function returning True or False (see below)
     * a collection of any combination of the above
     
    Return:
     * a list of filter functions (see below)
     
    Filter Functions:
     * for file searches, these take the path and filename as 2 separate arguments
     * for file contents, take a single line of text as one argument
    """
    
    if isinstance(filter_param, str) or isinstance(filter_param, re.Pattern):
        # return a list of just this pattern as a function
        return [_get_rexp_filter(filter_param, func_type),]
    
    elif callable(filter_param):
        # just put this function in a list
        return [filter_param,]
    
    else:
        # assume it's an iterable
        return [
            (
                fp 
                if callable(fp) 
                else _get_rexp_filter(fp, func_type)
            )
            for fp in filter_param
        ]


NOTEBOOK_FILTERS = _process_filter(
    [
        lambda path, name: _NB_REXP.search(name) is not None,
        lambda path, name: _CP_REXP.search(path) is None,
    ]
)

### Primary class

@attr.s
class NotebookSearcher():
    """
    Can search a tree of jupyter nobooks on both file names and file contents 

    Instantiate with your notebook dir:

        nb_searcher = NotebookSearcher('/path/to/notebooks')

    (The default is /home/user/notebooks).

    use nb_searcher.find_notebooks() to search.
    """

    notebook_root = attr.ib(default=f'{HOME}/notebooks')
    
    def find_notebooks(
        self,
        contents_filter=None,
        name_filter=None,
        keep_matches=False
    ):
        """
        returns notebooks that match filters

        params:
            {contents,name}_filter: single or list of strings or regular
            expressions to look for in either the file names or file contents
            keep_matches: if False, only return notebook names, otherwise keep
            detailed match locations
        """
        
        if name_filter is None:
            # default to non-checkpoint .ipynb files
            name_filter = NOTEBOOK_FILTERS
        else:
            # add more filename filters if passed
            name_filter = NOTEBOOK_FILTERS + _process_filter(name_filter)
        
        # option file contents filters
        if contents_filter is not None:
            contents_filter = _process_filter(contents_filter, func_type=SEARCH_CONTENTS)

        matched_notebooks = {}
        for notebook_file in filesystem.find(
            self.notebook_root, 
            filters=name_filter
        ):
            if contents_filter is None:
                matched_notebooks[notebook_file] = True
            else:
                with open(notebook_file) as lines:
                    matches = dict()
                    for i, line in enumerate(lines):
                        if keep_matches:
                            line_matches = []
                            for j, single_filter in enumerate(contents_filter):
                                if single_filter(line):
                                    line_matches.append(j)
                            if len(line_matches) > 0:
                                matches[i] = {'text': line, 'hits': line_matches}
                        else:
                            if any(single_filter(line) for single_filter in contents_filter):
                                matched_notebooks[notebook_file] = True
                                break
                    else:
                        if keep_matches and len(matches) > 0:
                            matched_notebooks[notebook_file] = matches
                    
        return matched_notebooks if keep_matches else set(matched_notebooks)

### Helper fns

def _get_rexp_filter(rexp, func_type=SEARCH_FILE):
    """
    Turn a pattern (string or re.Pattern) into a callable
    """
    if callable(rexp):
        if func_type == SEARCH_FILE:
            try:
                answer = rexp("path", "file")
                return rexp
            except TypeError:
                try:
                    answer = rexp("path/file")
                    return lambda p,f: rexp(os.path.join(p, f))
                except TypeError:
                    raise Exception("Function has to be able to take a string or path/filname pair!")
        elif func_type == SEARCH_CONTENTS:
            try:
                answer = rexp("lines")
                return rexp
            except TypeError:
                raise Exception("Function has to be able to take a single string ")
        else:
            raise Exception(f"I don't know what to do with search type: '{func_type}'")
            
    
    if isinstance(rexp, str):
        rexp = re.compile(rexp)
    elif not isinstance(rexp, re.Pattern):
        raise Exception("Filter must be a pattern (str or re.Pattern) or callable")
    
    if func_type == SEARCH_FILE:
        return lambda path, name: rexp.search(os.path.join(path, name)) is not None
    elif func_type == SEARCH_CONTENTS:
        return lambda line: rexp.search(line) is not None
    else:
        raise Exception(f"I don't know what to do with search type: '{func_type}'")

