import json
import os
import re

import pytest

from jme.jupy_tools import nb_searcher
from jme.jupy_tools.nb_searcher import NotebookSearcher


def _notebook_text(*cell_sources):
    """ minimal ipynb JSON, one cell source per line so line numbers are predictable """
    return json.dumps({
        "cells": [{"cell_type": "code", "source": [src]} for src in cell_sources],
        "metadata": {}, "nbformat": 4, "nbformat_minor": 5,
    }, indent=1)


@pytest.fixture
def nb_root(tmp_path):
    """
    tmp_path/
        assembly.ipynb          (uses minimap2 and pandas)
        notes.txt               (not a notebook, mentions minimap2)
        project_a/
            mapping.ipynb       (uses minimap2)
            plots.ipynb         (uses matplotlib)
            .ipynb_checkpoints/
                mapping-checkpoint.ipynb   (should always be ignored)
        project_b/
            deep/
                taxonomy.ipynb  (uses pandas)
    """
    files = {
        'assembly.ipynb': _notebook_text("!minimap2 a.fa b.fa", "import pandas"),
        'notes.txt': "remember to run minimap2\n",
        'project_a/mapping.ipynb': _notebook_text("!minimap2 -x map-ont"),
        'project_a/plots.ipynb': _notebook_text("from matplotlib import pyplot"),
        'project_a/.ipynb_checkpoints/mapping-checkpoint.ipynb':
            _notebook_text("!minimap2 -x map-ont"),
        'project_b/deep/taxonomy.ipynb': _notebook_text("import pandas as pd"),
    }
    for rel_path, text in files.items():
        path = tmp_path / rel_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    return tmp_path


def _rel(paths, root):
    return {os.path.relpath(p, root) for p in paths}


ALL_NOTEBOOKS = {'assembly.ipynb', 'project_a/mapping.ipynb',
                 'project_a/plots.ipynb', 'project_b/deep/taxonomy.ipynb'}


def test_default_root():
    assert NotebookSearcher().notebook_root == \
        os.path.join(os.environ['HOME'], 'notebooks')


def test_find_all_notebooks(nb_root):
    found = NotebookSearcher(str(nb_root)).find_notebooks()
    assert isinstance(found, set)
    assert _rel(found, nb_root) == ALL_NOTEBOOKS


def test_find_notebooks_no_matches(tmp_path):
    assert NotebookSearcher(str(tmp_path)).find_notebooks() == set()


@pytest.mark.parametrize("name_filter, expected", [
    # (patterns are matched against the full path, so anchor them to avoid
    #  matching pytest's tmp dir names)
    ("/project_a/", {'project_a/mapping.ipynb', 'project_a/plots.ipynb'}),
    (re.compile(r"/(mapping|taxonomy)\.ipynb$"),
     {'project_a/mapping.ipynb', 'project_b/deep/taxonomy.ipynb'}),
    (lambda path, name: name.startswith('a'), {'assembly.ipynb'}),
    # lists are ANDed
    (["/project_", re.compile("/plot")], {'project_a/plots.ipynb'}),
    (["/project_", lambda p, n: p.endswith('/deep')], {'project_b/deep/taxonomy.ipynb'}),
])
def test_name_filter(nb_root, name_filter, expected):
    found = NotebookSearcher(str(nb_root)).find_notebooks(name_filter=name_filter)
    assert _rel(found, nb_root) == expected


def test_name_filter_single_arg_callable(nb_root):
    """ one argument functions are called with the full path """
    found = NotebookSearcher(str(nb_root)).find_notebooks(
        name_filter=lambda full_path: '/deep/' in full_path)
    assert _rel(found, nb_root) == {'project_b/deep/taxonomy.ipynb'}

    # also in a list
    found = NotebookSearcher(str(nb_root)).find_notebooks(
        name_filter=["/project_", lambda full_path: full_path.endswith('s.ipynb')])
    assert _rel(found, nb_root) == {'project_a/plots.ipynb'}


def test_name_filter_uses_real_files(nb_root):
    """
    filters aren't called with made up paths to work out their arguments
    (this one would raise FileNotFoundError)
    """
    sizes = sorted(os.path.getsize(nb_root / nb) for nb in ALL_NOTEBOOKS)
    cutoff = sizes[0]

    def big_notebooks(path, name):
        return os.path.getsize(os.path.join(path, name)) > cutoff

    found = NotebookSearcher(str(nb_root)).find_notebooks(
        name_filter=big_notebooks)
    expected = {nb for nb in ALL_NOTEBOOKS
                if os.path.getsize(nb_root / nb) > cutoff}
    assert 0 < len(expected) < len(ALL_NOTEBOOKS)
    assert _rel(found, nb_root) == expected


@pytest.mark.parametrize("contents_filter, expected", [
    ("minimap2", {'assembly.ipynb', 'project_a/mapping.ipynb'}),
    (re.compile(r"import\s+pandas"),
     {'assembly.ipynb', 'project_b/deep/taxonomy.ipynb'}),
    (lambda line: 'matplotlib' in line, {'project_a/plots.ipynb'}),
    # multiple content filters are ORed
    (["matplotlib", "map-ont"],
     {'project_a/plots.ipynb', 'project_a/mapping.ipynb'}),
    ("not in any notebook", set()),
])
def test_contents_filter(nb_root, contents_filter, expected):
    found = NotebookSearcher(str(nb_root)).find_notebooks(
        contents_filter=contents_filter)
    assert _rel(found, nb_root) == expected


def test_name_and_contents_filter(nb_root):
    found = NotebookSearcher(str(nb_root)).find_notebooks(
        contents_filter="pandas", name_filter="/project_")
    assert _rel(found, nb_root) == {'project_b/deep/taxonomy.ipynb'}


def test_keep_matches(nb_root):
    found = NotebookSearcher(str(nb_root)).find_notebooks(
        contents_filter=["minimap2", "pandas"], keep_matches=True)
    assert isinstance(found, dict)
    found = {os.path.relpath(k, nb_root): v for k, v in found.items()}
    assert set(found) == {'assembly.ipynb', 'project_a/mapping.ipynb',
                          'project_b/deep/taxonomy.ipynb'}

    # line number -> text and indices of the filters that hit
    assembly = found['assembly.ipynb']
    assert len(assembly) == 2
    hits = sorted((v['hits'], v['text'].strip()) for v in assembly.values())
    assert hits[0][0] == [0]
    assert 'minimap2' in hits[0][1]
    assert hits[1][0] == [1]
    assert 'pandas' in hits[1][1]

    lines = _notebook_text("!minimap2 a.fa b.fa", "import pandas").split("\n")
    for line_no, match in assembly.items():
        assert lines[line_no] == match['text'].rstrip("\n")


def test_keep_matches_without_contents_filter(nb_root):
    found = NotebookSearcher(str(nb_root)).find_notebooks(keep_matches=True)
    assert isinstance(found, dict)
    assert _rel(found, nb_root) == ALL_NOTEBOOKS
    assert all(v is True for v in found.values())


###
# helpers


def test_notebook_filters():
    assert all(f('/some/path', 'nb.ipynb') for f in nb_searcher.NOTEBOOK_FILTERS)
    assert not all(f('/some/path', 'nb.py') for f in nb_searcher.NOTEBOOK_FILTERS)
    assert not all(f('/some/.ipynb_checkpoints', 'nb-checkpoint.ipynb')
                   for f in nb_searcher.NOTEBOOK_FILTERS)


def test_get_rexp_filter_file():
    file_filter = nb_searcher._get_rexp_filter("a/b")
    assert file_filter("x/a", "b.ipynb")
    assert not file_filter("x/b", "a.ipynb")

    file_filter = nb_searcher._get_rexp_filter(re.compile("^/x"),
                                               nb_searcher.SEARCH_FILE)
    assert file_filter("/x/y", "z")
    assert not file_filter("/y/x", "z")


def test_get_rexp_filter_contents():
    line_filter = nb_searcher._get_rexp_filter("foo",
                                               nb_searcher.SEARCH_CONTENTS)
    assert line_filter("a foo b")
    assert not line_filter("a bar b")


def test_get_rexp_filter_callables():
    two_arg = lambda path, name: name == 'x'
    assert nb_searcher._get_rexp_filter(two_arg) is two_arg

    one_arg = nb_searcher._get_rexp_filter(lambda full: full == 'a/x')
    assert one_arg('a', 'x')
    assert not one_arg('b', 'x')

    line_fn = lambda line: 'x' in line
    assert nb_searcher._get_rexp_filter(line_fn,
                                        nb_searcher.SEARCH_CONTENTS) is line_fn


def test_get_rexp_filter_not_called():
    """ the number of arguments comes from the signature, not a test call """
    calls = []

    def two_args(path, name):
        calls.append((path, name))
        return True

    def one_arg(full_path):
        calls.append(full_path)
        return True

    assert nb_searcher._get_rexp_filter(two_args) is two_args
    wrapped = nb_searcher._get_rexp_filter(one_arg)
    assert nb_searcher._get_rexp_filter(one_arg,
                                        nb_searcher.SEARCH_CONTENTS) is one_arg
    assert calls == []
    assert wrapped('a', 'b')
    assert calls == ['a/b']


def test_get_rexp_filter_internal_type_error():
    """ a TypeError inside a 2 argument filter isn't mistaken for 1 argument """
    def bad_filter(path, name):
        return name + 1

    assert nb_searcher._get_rexp_filter(bad_filter) is bad_filter


def test_get_rexp_filter_other_callables():
    import functools

    def has_ext(path, name, ext):
        return name.endswith(ext)

    notebooks = functools.partial(has_ext, ext='.ipynb')
    assert nb_searcher._get_rexp_filter(notebooks) is notebooks

    class FullPathFilter:
        def __call__(self, full_path):
            return 'x' in full_path

    adapted = nb_searcher._get_rexp_filter(FullPathFilter())
    assert adapted('a', 'x') and not adapted('a', 'b')

    # optional and variable arguments
    flexible = lambda *args: True
    assert nb_searcher._get_rexp_filter(flexible) is flexible
    with_default = lambda full_path, extra=None: True
    assert nb_searcher._get_rexp_filter(with_default) is with_default


class NoSignature:
    """ a callable whose signature can't be read (like some builtins) """
    __signature__ = "not a signature"

    def __init__(self, func):
        self.func = func

    def __call__(self, *args):
        return self.func(*args)


def test_get_rexp_filter_no_signature():
    """ without a signature, fall back to a test call """
    two_args = NoSignature(lambda path, name: name == 'x')
    assert nb_searcher._get_rexp_filter(two_args) is two_args

    one_arg = nb_searcher._get_rexp_filter(NoSignature(lambda full: full == 'a/x'))
    assert one_arg('a', 'x') and not one_arg('b', 'x')

    line_fn = NoSignature(lambda line: 'x' in line)
    assert nb_searcher._get_rexp_filter(line_fn,
                                        nb_searcher.SEARCH_CONTENTS) is line_fn
    with pytest.raises(Exception, match="path/filname pair"):
        nb_searcher._get_rexp_filter(NoSignature(lambda: True))
    with pytest.raises(Exception, match="single string"):
        nb_searcher._get_rexp_filter(NoSignature(lambda a, b: True),
                                     nb_searcher.SEARCH_CONTENTS)


def test_get_rexp_filter_errors():
    with pytest.raises(Exception, match="must be a pattern"):
        nb_searcher._get_rexp_filter(42)
    with pytest.raises(Exception, match="search type"):
        nb_searcher._get_rexp_filter("x", "bad_type")
    with pytest.raises(Exception, match="search type"):
        nb_searcher._get_rexp_filter(lambda x: x, "bad_type")
    with pytest.raises(Exception, match="path/filname pair"):
        nb_searcher._get_rexp_filter(lambda: True)
    with pytest.raises(Exception, match="single string"):
        nb_searcher._get_rexp_filter(lambda a, b: True,
                                     nb_searcher.SEARCH_CONTENTS)


def test_process_filter():
    fn = lambda p, n: True
    assert nb_searcher._process_filter(fn) == [fn]

    filters = nb_searcher._process_filter("abc")
    assert len(filters) == 1
    assert filters[0]("x", "abc")

    filters = nb_searcher._process_filter(["a", re.compile("b"), fn])
    assert len(filters) == 3
    assert filters[2] is fn
    assert filters[0]("a", "") and not filters[0]("b", "")
    assert filters[1]("b", "") and not filters[1]("a", "")

    line_filters = nb_searcher._process_filter(["a"], nb_searcher.SEARCH_CONTENTS)
    assert line_filters[0]("cat")
