from jme.jupy_tools import filesystem

def test_find():
    test_scripts = list(filesystem.find('tests', [filesystem.get_rexp_filter(r'\.py$')])) 
    assert 'tests/fs_test.py' in test_scripts
    assert 'tests/utils_test.py' in test_scripts
    assert 'tests/cdhit_test.py' in test_scripts
    assert 'tests/test.yml' not in test_scripts

def test_glob_wildcards():
    results = {f:w for f,w in filesystem.glob_wildcards('tests/{mod}_test.py')}
    assert 'tests/fs_test.py' in results
    assert 'tests/utils_test.py' in results
    assert 'tests/cdhit_test.py' in results
    assert 'tests/test.yml' not in results

    assert results['tests/fs_test.py']['mod'] == 'fs'
    

def test_glob_wildcards_tuple():
    results = {f:w for f,w in filesystem.glob_wildcards('tests/{mod}_{suff}.py', as_tuple=True)}
    assert 'tests/fs_test.py' in results
    assert 'tests/utils_test.py' in results
    assert 'tests/cdhit_test.py' in results
    assert 'tests/test.yml' not in results
    fs_wc = results['tests/fs_test.py']
    mod, suff = fs_wc
    assert mod == fs_wc.mod
    assert suff == fs_wc.suff
    assert suff == "test"

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
