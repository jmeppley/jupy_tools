import importlib
import sys
import warnings

import pytest


def _star_import():
    """ run `from jme.jupy_tools.faves import *` in a fresh namespace """
    sys.modules.pop('jme.jupy_tools.faves', None)
    namespace = {}
    exec("from jme.jupy_tools.faves import *", namespace)
    return namespace


@pytest.fixture(autouse=True)
def _restore_faves():
    yield
    sys.modules.pop('jme.jupy_tools.faves', None)


def test_star_import():
    from jme.jupy_tools import conda, hit_tables, filesystem, utils
    import collections, itertools, os, re, json
    import numpy, pandas
    from matplotlib import pyplot

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ns = _star_import()

    # from this package
    assert ns['conda'] is conda
    assert ns['hit_tables'] is hit_tables
    assert ns['glob_wildcards'] is filesystem.glob_wildcards
    assert ns['go_to_dir'] is filesystem.go_to_dir
    assert ns['go2dir'] is filesystem.go_to_dir
    for name in ['LogLogger', 'first', 'read_tsv', 'save_fig_to_pdf']:
        assert ns[name] is getattr(utils, name)

    # standard library
    assert ns['Counter'] is collections.Counter
    assert ns['defaultdict'] is collections.defaultdict
    for name in ['islice', 'combinations', 'chain']:
        assert ns[name] is getattr(itertools, name)
    for module in [os, re, json, sys, warnings]:
        assert ns[module.__name__] is module

    # required third party
    assert ns['np'] is numpy
    assert ns['pd'] is pandas
    assert ns['plt'] is pyplot


def test_optional_imports():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ns = _star_import()
    for name, module in [('yaml', 'yaml'), ('SeqIO', 'Bio'),
                         ('sns', 'seaborn'), ('pl', 'polars')]:
        try:
            importlib.import_module(module)
        except ImportError:
            assert name not in ns
        else:
            assert name in ns
    if 'SeqIO' in ns:
        from jme.jupy_tools.biopython import get_seq_lens
        assert ns['get_seq_lens'] is get_seq_lens


@pytest.mark.parametrize("module, message", [
    ('yaml', 'yaml'),
    ('Bio', 'biopython'),
    ('seaborn', 'seaborn'),
    ('polars', 'Polars'),
])
def test_missing_optional_module_warns(monkeypatch, module, message):
    monkeypatch.setitem(sys.modules, module, None)
    with pytest.warns(UserWarning, match=f"{message} import failed"):
        _star_import()
