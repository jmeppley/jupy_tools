import os
import time
import warnings
from datetime import datetime

import numpy
import pandas
import pytest

from jme.jupy_tools import filesystem


@pytest.fixture
def file_tree(tmp_path):
    """
    tmp_path/
        a.txt
        b.py
        sub/
            c.txt
            d.py
            deeper/
                e.txt
        skipme/
            f.txt
    """
    for rel_path in ['a.txt', 'b.py', 'sub/c.txt', 'sub/d.py',
                     'sub/deeper/e.txt', 'skipme/f.txt']:
        path = tmp_path / rel_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(rel_path)
    return tmp_path


###
# find


def test_find():
    test_scripts = list(filesystem.find('tests', [filesystem.get_rexp_filter(r'\.py$')]))
    assert 'tests/test_fs.py' in test_scripts
    assert 'tests/test_utils.py' in test_scripts
    assert 'tests/test_cdhit.py' in test_scripts
    assert 'tests/test.yml' not in test_scripts


def test_find_no_filters(file_tree):
    found = {os.path.relpath(p, file_tree) for p in filesystem.find(str(file_tree))}
    assert found == {'a.txt', 'b.py', 'sub', 'sub/c.txt', 'sub/d.py',
                     'sub/deeper', 'sub/deeper/e.txt', 'skipme', 'skipme/f.txt'}


def test_find_filters_are_anded(file_tree):
    filters = [filesystem.get_rexp_filter(r'\.txt$'),
               filesystem.get_rexp_filter('sub')]
    found = {os.path.relpath(p, file_tree)
             for p in filesystem.find(str(file_tree), filters)}
    assert found == {'sub/c.txt', 'sub/deeper/e.txt'}


def test_find_custom_filter(file_tree):
    def no_skipme(root, name):
        return 'skipme' not in os.path.join(root, name)
    found = {os.path.relpath(p, file_tree)
             for p in filesystem.find(str(file_tree), [no_skipme])}
    assert 'skipme' not in found
    assert 'skipme/f.txt' not in found
    assert 'sub/deeper/e.txt' in found


def test_get_rexp_filter_compiled():
    import re
    rexp_filter = filesystem.get_rexp_filter(re.compile(r'^/data/.+\.fa$'))
    assert rexp_filter('/data/genomes', 'x.fa')
    assert not rexp_filter('/data/genomes', 'x.fa.gz')
    assert not rexp_filter('/other', 'x.fa')


def test_check_path():
    assert filesystem._check_path('root', 'name')
    assert filesystem._check_path('root', 'name', [])
    assert filesystem._check_path('root', 'name', [lambda r, n: True])
    assert not filesystem._check_path('root', 'name',
                                      [lambda r, n: True, lambda r, n: False])


###
# glob_wildcards


def test_glob_wildcards():
    results = {f: w for f, w in filesystem.glob_wildcards('tests/test_{mod}.py')}
    assert 'tests/test_fs.py' in results
    assert 'tests/test_utils.py' in results
    assert 'tests/test_cdhit.py' in results
    assert 'tests/test.yml' not in results

    assert results['tests/test_fs.py']['mod'] == 'fs'


def test_glob_wildcards_tuple():
    results = {f: w for f, w in filesystem.glob_wildcards('tests/{pref}_{mod}.py', return_type=tuple)}
    assert 'tests/test_fs.py' in results
    assert 'tests/test_utils.py' in results
    assert 'tests/test_cdhit.py' in results
    assert 'tests/test.yml' not in results
    fs_wc = results['tests/test_fs.py']
    pref, mod = fs_wc
    assert mod == fs_wc.mod
    assert pref == fs_wc.pref
    assert pref == "test"
    assert mod == "fs"


@pytest.fixture
def sample_files(tmp_path):
    names = ['s1.reads.fastq', 's2.reads.fastq', 's2.reads.fasta',
             's3.contigs.fasta', 's10.reads.fastq', 'other.txt', 's1.s1.txt',
             's1.s2.txt']
    for name in names:
        (tmp_path / name).write_text(name)
    return tmp_path


def test_glob_wildcards_constraints(sample_files):
    template = str(sample_files / '{sample}.reads.{ext}')
    all_hits = dict(filesystem.glob_wildcards(template))
    assert len(all_hits) == 4

    hits = dict(filesystem.glob_wildcards(template,
                                          constraints={'sample': r's\d'}))
    assert {w['sample'] for w in hits.values()} == {'s1', 's2'}
    assert {w['ext'] for w in hits.values()} == {'fastq', 'fasta'}


def test_glob_wildcards_dots_are_literal(sample_files):
    # the "." in the template must not match any character
    template = str(sample_files / '{sample}.reads.fastq')
    (sample_files / 's4Xreads.fastq').write_text('nope')
    samples = {w['sample'] for _, w in filesystem.glob_wildcards(template)}
    assert samples == {'s1', 's2', 's10'}


def test_glob_wildcards_repeated_name(sample_files):
    """ a wildcard used twice must match the same value both times """
    template = str(sample_files / '{sample}.{sample}.txt')
    hits = dict(filesystem.glob_wildcards(template))
    assert list(hits) == [str(sample_files / 's1.s1.txt')]
    assert hits[str(sample_files / 's1.s1.txt')] == {'sample': 's1'}


def test_glob_wildcards_snek(sample_files):
    template = str(sample_files / '{sample}.reads.{ext}')
    samples, exts = filesystem.glob_wildcards(template,
                                              return_type=filesystem.SNEK)
    assert sorted(samples) == ['s1', 's10', 's2', 's2']
    assert sorted(exts) == ['fasta', 'fastq', 'fastq', 'fastq']
    # values stay paired up
    assert sorted(zip(samples, exts)) == [('s1', 'fastq'), ('s10', 'fastq'),
                                          ('s2', 'fasta'), ('s2', 'fastq')]


@pytest.mark.parametrize('return_type', [filesystem.TUPL, tuple,
                                         __import__('collections').namedtuple])
def test_glob_wildcards_tuple_aliases(sample_files, return_type):
    template = str(sample_files / '{sample}.contigs.{ext}')
    (path, wildcards), = filesystem.glob_wildcards(template,
                                                    return_type=return_type)
    assert path == str(sample_files / 's3.contigs.fasta')
    assert wildcards.sample == 's3'
    assert wildcards.ext == 'fasta'


def test_glob_wildcards_dict_alias(sample_files):
    template = str(sample_files / '{sample}.contigs.fasta')
    (path, wildcards), = filesystem.glob_wildcards(template, return_type=dict)
    assert wildcards == {'sample': 's3'}


def test_glob_wildcards_as_tuple_deprecated(sample_files):
    template = str(sample_files / '{sample}.contigs.fasta')
    with pytest.warns(UserWarning, match='return_type'):
        (path, wildcards), = filesystem.glob_wildcards(template, as_tuple=True)
    assert wildcards.sample == 's3'

    # as_tuple is ignored if return_type is given
    with pytest.warns(UserWarning, match='ignored'):
        samples, = filesystem.glob_wildcards(template, as_tuple=True,
                                             return_type=filesystem.SNEK)
    assert samples == ['s3']


def test_glob_wildcards_no_matches(tmp_path):
    assert list(filesystem.glob_wildcards(str(tmp_path / '{x}.nothing'))) == []
    # SNEK mode yields one empty list per wildcard
    assert list(filesystem.glob_wildcards(str(tmp_path / '{x}.{y}.nothing'),
                                          return_type=filesystem.SNEK)) == [[], []]


def test_glob_wildcards_debug(sample_files, capsys):
    template = str(sample_files / '{sample}.reads.fastq')
    list(filesystem.glob_wildcards(template, debug=True,
                                   constraints={'sample': r's\d'}))
    out = capsys.readouterr().out
    assert 'Wildcard regexp' in out
    # s10 is found by glob, but rejected by the constraint
    assert "s10.reads.fastq doesn match" in out


def test_template_rexp_ignores_double_brackets():
    assert filesystem.TEMPLATE_REXP.findall("{a}/{{b}}/{ c }") == ['a', 'c']


def test_hide_dots():
    assert filesystem._hide_dots("a.b.c") == r"a\.b\.c"


###
# human_readable_bytes


def test_human_readable_bytes():
    for val, code in [
        (2, "2"),
        (54, "54"),
        (243, ".24K"),
        (59049, "58K"),
        (4782969, "4.6M"),
        (7625597484987, "6.9T"),
        (150094635296999121, ".13E")
    ]:
        assert code == filesystem.human_readable_bytes(val)


@pytest.mark.parametrize("val, code", [
    (0, "0"),
    (99, "99"),
    (1024, "1.0K"),
    (1024 ** 2, "1.0M"),
    (5 * 1024 ** 3, "5.0G"),
    (1024 ** 5, "1.0P"),
])
def test_human_readable_bytes_more(val, code):
    assert filesystem.human_readable_bytes(val) == code


###
# usage tables


def test_get_user_lookup_table():
    users = filesystem.get_user_lookup_table()
    assert isinstance(users, pandas.Series)
    # root is always uid 0
    assert users[0] == 'root'


def test_get_bin_bounds_string():
    bounds = [0, 5, 10, 15]
    assert filesystem.get_bin_bounds_string(1, bounds) == "5 to 10 "
    assert filesystem.get_bin_bounds_string(2, bounds, str, "months") == \
        "10 to 15 months"
    assert filesystem.get_bin_bounds_string(0, bounds, lambda b: f"<{b}>",
                                            "x") == "<0> to <5> x"


def _make_aged_files(root):
    """ create 3 files: one new, one 10 days old, one 20 days old """
    now = time.time()
    day = 3600 * 24
    sizes = {}
    for name, age_days, size in [('new.txt', 0, 10),
                                 ('ten.txt', 10, 100),
                                 ('twenty.txt', 20, 1000)]:
        path = root / name
        path.write_bytes(b"x" * size)
        mtime = now - age_days * day
        os.utime(path, (mtime, mtime))
        sizes[name] = size
    return now, sizes


def test_get_file_sizes_and_dates_by_uid(tmp_path):
    now, sizes = _make_aged_files(tmp_path)
    # a broken link should be skipped
    os.symlink(tmp_path / 'missing', tmp_path / 'broken_link')

    usage, min_date = filesystem.get_file_sizes_and_dates_by_uid(str(tmp_path))
    assert len(usage) == 1
    owner, = usage
    entries = usage[owner]
    assert sorted(e[0] for e in entries) == [10, 100, 1000]
    assert {os.path.basename(e[2]) for e in entries} == set(sizes)
    assert abs(min_date - (now - 20 * 3600 * 24)) < 5

    # min_age drops new files (but min_date still reflects all files)
    usage, min_date_2 = filesystem.get_file_sizes_and_dates_by_uid(
        str(tmp_path), min_age=5 * 3600 * 24)
    assert sorted(e[0] for e in usage[owner]) == [100, 1000]
    assert min_date_2 == min_date


def test_get_file_sizes_and_dates_by_uid_users(tmp_path):
    _make_aged_files(tmp_path)
    me = os.getuid()
    usage, _ = filesystem.get_file_sizes_and_dates_by_uid(str(tmp_path),
                                                          users=[me])
    assert sum(len(v) for v in usage.values()) == 3

    usage, _ = filesystem.get_file_sizes_and_dates_by_uid(
        str(tmp_path), users=['no_such_user_hopefully'])
    assert len(usage) == 0


@pytest.mark.parametrize('time_stat', ['mtime', 'ctime', 'atime', 'max'])
def test_get_file_sizes_and_dates_by_uid_time_stat(tmp_path, time_stat):
    _make_aged_files(tmp_path)
    usage, min_date = filesystem.get_file_sizes_and_dates_by_uid(
        str(tmp_path), time_stat=time_stat)
    assert sum(len(v) for v in usage.values()) == 3


def test_get_file_size_table(tmp_path):
    _make_aged_files(tmp_path)
    usage, min_date = filesystem.get_file_sizes_and_dates_by_uid(str(tmp_path))
    owner, = usage

    table = filesystem.get_file_size_table(usage, min_date, age_bin_size=1,
                                           age_bin_type='weeks')
    assert list(table.columns) == [owner]
    assert table[owner].sum() == 1110
    sizes_by_bin = table[owner].to_dict()
    assert sizes_by_bin['0 to 1 weeks old'] == 10
    assert sizes_by_bin['1 to 2 weeks old'] == 100
    assert sizes_by_bin['2 to 3 weeks old'] == 1000

    # min_age skips the newest file
    table = filesystem.get_file_size_table(usage, min_date, age_bin_size=1,
                                           age_bin_type='weeks',
                                           min_age=3600 * 24)
    assert table[owner].sum() == 1100


def test_get_file_size_table_bad_span():
    with pytest.raises(Exception, match="fortnights"):
        filesystem.get_file_size_table({}, 0, age_bin_type='fortnights')


def test_time_spans():
    assert filesystem.TIME_SPANS['days'] == 24 * filesystem.TIME_SPANS['hours']
    assert filesystem.TIME_SPANS['weeks'] == 7 * filesystem.TIME_SPANS['days']


###
# go_to_dir


def test_go_to_dir(tmp_path, capsys):
    new_dir = tmp_path / 'some' / 'new' / 'dir'
    assert filesystem.go_to_dir(str(new_dir)) == str(new_dir)
    assert os.path.isdir(new_dir)
    assert os.path.samefile(os.getcwd(), new_dir)
    assert capsys.readouterr().out.strip() == str(new_dir)

    # existing dir, quietly
    assert filesystem.go_to_dir(str(tmp_path), quiet=True) == str(tmp_path)
    assert os.path.samefile(os.getcwd(), tmp_path)
    assert capsys.readouterr().out == ''
