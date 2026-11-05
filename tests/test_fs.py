from jme.jupy_tools import filesystem

def test_find():
    test_scripts = list(filesystem.find('tests', [filesystem.get_rexp_filter(r'\.py$')])) 
    assert 'tests/test_fs.py' in test_scripts
    assert 'tests/test_utils.py' in test_scripts
    assert 'tests/test_cdhit.py' in test_scripts
    assert 'tests/test.yml' not in test_scripts

def test_glob_wildcards():
    results = {f:w for f,w in filesystem.glob_wildcards('tests/test_{mod}.py')}
    assert 'tests/test_fs.py' in results
    assert 'tests/test_utils.py' in results
    assert 'tests/test_cdhit.py' in results
    assert 'tests/test.yml' not in results

    assert results['tests/test_fs.py']['mod'] == 'fs'
    

def test_glob_wildcards_tuple():
    results = {f:w for f,w in filesystem.glob_wildcards('tests/{pref}_{mod}.py', return_type=tuple)}
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
