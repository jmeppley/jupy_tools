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
    

