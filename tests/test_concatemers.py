import sys

import numpy
import pandas
import pytest

from jme.jupy_tools import concatemers, hit_tables


def tandem_self_hits(read, length, repeat_size, pctid=99.0):
    """
    Build the self-vs-self hit table (BLAST columns) for a perfect tandem
    repeat: one hit for each offset k * repeat_size, including the full length
    self hit at k=0.
    """
    rows = []
    k_max = (length - 1) // repeat_size
    for k in range(-k_max, k_max + 1):
        mlen = length - abs(k) * repeat_size
        if k >= 0:
            qstart, qend, hstart, hend = 1, mlen, 1 + k * repeat_size, length
        else:
            qstart, qend, hstart, hend = 1 - k * repeat_size, length, 1, mlen
        rows.append([read, read, pctid, mlen, 0, 0, qstart, qend, hstart, hend,
                     0.0, float(mlen)])
    return pandas.DataFrame(rows, columns=hit_tables.BLAST_COLUMNS)


def to_m8(hits):
    return hits.to_csv(sep="\t", header=False, index=False)


@pytest.fixture
def repeat_hits():
    """ a 10kb read made of 5 copies of a 2kb repeat """
    return tandem_self_hits('r1', 10000, 2000)


@pytest.fixture
def single_hits():
    """ a read with only its full length self hit """
    return tandem_self_hits('r2', 5000, 5000)


@pytest.fixture
def scipy():
    return pytest.importorskip("scipy")


@pytest.fixture
def sklearn():
    return pytest.importorskip("sklearn")


###
# basic calculations


def test_get_ranks(repeat_hits):
    ranks = concatemers.get_ranks(repeat_hits)
    assert sorted(ranks) == [k * 2000 for k in range(-4, 5)]


def test_get_rank_sizes(repeat_hits):
    rank_sizes = concatemers.get_rank_sizes(repeat_hits)
    assert rank_sizes[0] == 10000
    assert rank_sizes[2000] == 8000
    assert rank_sizes[-8000] == 2000
    assert len(rank_sizes) == 9


def test_get_rank_sizes_sums_fragments():
    hits = pandas.DataFrame({
        'qstart': [1, 501], 'qend': [100, 600],
        'hstart': [11, 511], 'hend': [110, 610],
        'mlen': [100, 100],
    })
    assert dict(concatemers.get_rank_sizes(hits)) == {10: 200}


def test_weighted_ranks():
    weighted = concatemers.weighted_ranks({0: 3000, 100: 1500, 200: 700})
    assert sorted(weighted) == [0, 0, 0, 0, 100, 100]
    assert len(concatemers.weighted_ranks({0: 3000}, factor=1000)) == 3


def test_evaluate_repeat():
    assert concatemers.evaluate_repeat(2000, 0.9, 10000)
    # score too low
    assert not concatemers.evaluate_repeat(2000, 0.4, 10000)
    # repeat is too big a fraction of the read
    assert not concatemers.evaluate_repeat(8000, 0.9, 10000)
    # custom cutoffs
    assert concatemers.evaluate_repeat(8000, 0.4, 10000, score_cutoff=0.3,
                                       frac_cutoff=0.8)


def test_get_seq_length(repeat_hits):
    assert concatemers.get_seq_length(repeat_hits, 1234) == 1234
    # from the largest qend if there is no qlen column
    assert concatemers.get_seq_length(repeat_hits, None) == 10000
    with_qlen = repeat_hits.assign(qlen=12000)
    assert concatemers.get_seq_length(with_qlen, None) == 12000


@pytest.mark.parametrize("hit_range, size, expected", [
    ((3, 5), 8, [0, 0, 1, 1, 1, 0, 0, 0]),
    ((5, 3), 8, [0, 0, 1, 1, 1, 0, 0, 0]),
    ((0, 2), 5, [1, 1, 0, 0, 0]),
    ((4, 99), 5, [0, 0, 0, 1, 1]),
    ((1, 5), 5, [1, 1, 1, 1, 1]),
    # floats are truncated
    ((2.7, 3.9), 5, [0, 1, 1, 0, 0]),
])
def test_make_presence_array(hit_range, size, expected):
    assert list(concatemers.make_presence_array(hit_range, size)) == expected


def test_get_score(repeat_hits):
    assert concatemers.get_score(repeat_hits, 2000) == pytest.approx(1.0)
    assert concatemers.get_score(repeat_hits, 2000, qlen=10000) == \
        pytest.approx(1.0)


def test_get_score_bad_guess(repeat_hits):
    # no hits near this offset
    assert concatemers.get_score(repeat_hits, 3000) == 0
    # only covers half the expected range
    half = repeat_hits[repeat_hits.qend != 8000]
    half = pandas.concat([half, pandas.DataFrame([{
        'qstart': 1, 'qend': 4000, 'hstart': 2001, 'hend': 6000,
        'mlen': 4000}])])
    assert concatemers.get_score(half, 2000, qlen=10000) == pytest.approx(0.5)


def test_generate_xys_from_clusters():
    xys = list(concatemers.generate_xys_from_clusters({0: [1, 2], 5: [3]}))
    assert xys == [(0, 1), (0, 2), (5, 3)]


def test_check_repeat_size(repeat_hits):
    ranks = concatemers.get_ranks(repeat_hits).values
    found, expected = concatemers.check_repeat_size(ranks, 2000, 10000,
                                                    refine=False)
    assert found == 9

    found, expected, (slope, intercept) = \
        concatemers.check_repeat_size(ranks, 2000, 10000)
    assert found == 9
    assert slope == pytest.approx(2000)
    assert intercept == 0

    # refinement corrects a slightly wrong guess
    found, expected, (slope, intercept) = \
        concatemers.check_repeat_size(ranks, 1990, 10000)
    assert slope == pytest.approx(2000)


def test_check_repeat_size_epsilon(repeat_hits):
    ranks = concatemers.get_ranks(repeat_hits).values
    # only ranks within epsilon/2 of a multiple are kept (plus zero)
    found, _ = concatemers.check_repeat_size(ranks, 1900, 10000, refine=False,
                                             epsilon=100)
    assert found == 1


@pytest.mark.xfail(strict=True, reason=(
    "check_repeat_size computes the expected cluster count from "
    "repeat_size / qlen instead of qlen / repeat_size"))
def test_check_repeat_size_expected_clusters(repeat_hits):
    ranks = concatemers.get_ranks(repeat_hits).values
    found, expected = concatemers.check_repeat_size(ranks, 2000, 10000,
                                                    refine=False)
    assert expected == found


###
# repeat size finders


def test_find_cmer_stats_fft(scipy, repeat_hits):
    repeat_size = concatemers.find_cmer_stats_fft(repeat_hits)
    assert repeat_size == pytest.approx(2000, rel=0.01)

    repeat_size = concatemers.find_cmer_stats_fft(
        tandem_self_hits('r', 7000, 3000))
    assert repeat_size == pytest.approx(3000, rel=0.01)


def test_find_cmer_stats_fft_debug(scipy, repeat_hits):
    repeat_size, rank_sizes, X, signal, fft, freq = \
        concatemers.find_cmer_stats_fft(repeat_hits, debug=True,
                                        resolution=1000)
    assert repeat_size == pytest.approx(2000, rel=0.01)
    assert len(rank_sizes) == 9
    assert len(X) == 1000
    assert len(fft) == len(freq) == 500


def test_find_cmer_stats_fft_no_repeat(scipy, single_hits):
    assert concatemers.find_cmer_stats_fft(single_hits) is None


def test_find_cmer_stats_fft_ignores_reverse_hits(scipy, repeat_hits):
    reverse = pandas.DataFrame([{
        'query': 'r1', 'hit': 'r1', 'pctid': 99, 'mlen': 5000,
        'qstart': 1, 'qend': 5000, 'hstart': 9000, 'hend': 4001}])
    hits = pandas.concat([repeat_hits, reverse])
    assert concatemers.find_cmer_stats_fft(hits) == \
        concatemers.find_cmer_stats_fft(repeat_hits)


def test_find_cmer_stats_clust(sklearn, repeat_hits):
    assert concatemers.find_cmer_stats_clust(repeat_hits) == \
        pytest.approx(2000)


def test_find_cmer_stats_clust_no_repeat(sklearn, single_hits):
    assert concatemers.find_cmer_stats_clust(single_hits) is None


def test_find_cmer_stats_clust_too_few_clusters(sklearn):
    # two copies: only 3 ranks, which isn't enough after the hit_cutoff
    hits = tandem_self_hits('r', 4000, 2000)
    assert concatemers.find_cmer_stats_clust(hits, hit_cutoff=500) is None


def test_find_cmer_stats_jb_too_few_hits(single_hits):
    assert concatemers.find_cmer_stats_jb(single_hits) is None


@pytest.mark.xfail(strict=True, raises=NameError, reason=(
    "find_cmer_stats_jb uses np.median, but numpy is imported as 'numpy'"))
def test_find_cmer_stats_jb(repeat_hits):
    repeat_size = concatemers.find_cmer_stats_jb(repeat_hits)
    assert repeat_size is not None


###
# check_cmer


def test_check_cmer_fft(scipy, repeat_hits):
    repeat_size, copies, score, state, qlen = \
        concatemers.check_cmer(repeat_hits)
    assert repeat_size == pytest.approx(2000, rel=0.01)
    assert copies == pytest.approx(5, rel=0.01)
    assert score == pytest.approx(1, rel=0.01)
    assert state
    assert qlen == 10000


def test_check_cmer_clust(sklearn, repeat_hits):
    result = concatemers.check_cmer(repeat_hits, method=concatemers.CLUST)
    assert result == (2000, 5, 1, True, 10000)
    # method names are case insensitive
    assert concatemers.check_cmer(repeat_hits, method='CLUSTER') == result


def test_check_cmer_no_repeat(scipy, single_hits):
    assert concatemers.check_cmer(single_hits) == (None, None, -1, False, 5000)


def test_check_cmer_bad_method(repeat_hits):
    with pytest.raises(Exception, match="Unkown cmer method"):
        concatemers.check_cmer(repeat_hits, method='magic')


@pytest.mark.xfail(strict=True, reason=(
    "check_cmer calls find_cmer_stats_clust when method is JB"))
def test_check_cmer_jb(repeat_hits, monkeypatch):
    monkeypatch.setattr(concatemers, 'find_cmer_stats_jb',
                        lambda qhits, **kw: 2000)
    monkeypatch.setattr(concatemers, 'find_cmer_stats_clust',
                        lambda qhits, **kw: 1234)
    assert concatemers.check_cmer(repeat_hits, method=concatemers.JB)[0] == 2000


def test_check_cmer_from_file(scipy, repeat_hits, tmp_path):
    path = tmp_path / "self.m8"
    path.write_text(to_m8(repeat_hits))
    repeat_size, copies, score, state, qlen = concatemers.check_cmer(str(path))
    assert repeat_size == pytest.approx(2000, rel=0.01)
    assert state


###
# get_qhits


def test_get_qhits(tmp_path, repeat_hits):
    reverse = pandas.DataFrame([{
        'query': 'r1', 'hit': 'r1', 'pctid': 99, 'mlen': 5000, 'mismatches': 0,
        'gaps': 0, 'qstart': 1, 'qend': 5000, 'hstart': 9000, 'hend': 4001,
        'evalue': 0, 'score': 1}])
    path = tmp_path / "self.m8"
    path.write_text(to_m8(pandas.concat([repeat_hits, reverse])))

    qhits = concatemers.get_qhits(str(path))
    # reverse hit is dropped
    assert qhits.shape[0] == 9

    # cutoffs
    assert concatemers.get_qhits(str(path), mlen_cutoff=6000).shape[0] == 5
    assert concatemers.get_qhits(str(path), pctid_cutoff=99.5).shape[0] == 0


def test_get_qhits_blast_plus(tmp_path, repeat_hits, capsys):
    hits = repeat_hits.assign(qlen=10000, hlen=10000, raw_score=1)
    path = tmp_path / "self.m8"
    path.write_text(to_m8(hits))

    qhits = concatemers.get_qhits(str(path), table_format=hit_tables.BLAST_PLUS)
    assert 'qlen' in qhits.columns

    # deprecated option
    qhits = concatemers.get_qhits(str(path), blast_plus=True)
    assert 'qlen' in qhits.columns
    assert "DEPRECATED" in capsys.readouterr().out


def test_get_qhits_blast_plus_conflict(tmp_path, capsys):
    path = tmp_path / "self.paf"
    path.write_text(
        "r1\t1000\t0\t1000\t+\tr1\t1000\t0\t1000\t1000\t1000\t60\n"
        "r1\t1000\t0\t500\t-\tr1\t1000\t500\t1000\t500\t500\t60\n"
    )
    qhits = concatemers.get_qhits(str(path), table_format=hit_tables.PAF,
                                  blast_plus=True)
    assert "WARNING: ignoring blast_plus" in capsys.readouterr().out
    # PAF tables are filtered by strand
    assert list(qhits.strand) == ['+']


def test_insert_self_hit_bad_format(repeat_hits):
    with pytest.raises(Exception, match="self hit"):
        concatemers.insert_self_hit(repeat_hits, format=hit_tables.BLAST)


@pytest.mark.xfail(strict=True, raises=AttributeError, reason=(
    "insert_self_hit uses DataFrame.append, which was removed in pandas 2"))
def test_insert_self_hit():
    qhits = pandas.DataFrame([{
        'query': 'r1', 'qlen': 1000, 'qstart': 1, 'qend': 500, 'strand': '+',
        'hit': 'r1', 'hlen': 1000, 'hstart': 501, 'hend': 1000,
        'matches': 500, 'mlen': 500, 'mapqv': 0, 'pctid': 100}])
    with_self = concatemers.insert_self_hit(qhits)
    assert with_self.shape[0] == 2
    assert with_self.loc[-1, 'qend'] == 1000


###
# tables and scripting


@pytest.fixture
def all_v_all(repeat_hits, single_hits):
    other = tandem_self_hits('r3', 7000, 3000)
    # a hit between two different reads
    cross = pandas.DataFrame([[
        'r1', 'r3', 90.0, 1000, 0, 0, 1, 1000, 1, 1000, 0.0, 1000.0]],
        columns=hit_tables.BLAST_COLUMNS)
    return pandas.concat([repeat_hits, single_hits, other, cross],
                         ignore_index=True)


def test_build_all_v_all_cmer_table(scipy, all_v_all, tmp_path, capsys):
    out_file = tmp_path / "cmers.tsv"
    table = concatemers.build_all_v_all_cmer_table(all_v_all, str(out_file))
    assert list(table.columns) == ['read', 'repeat_size', 'copies',
                                   'repeat_score', 'state', 'qlen']
    # r2 has no repeat and is dropped by the default 'notna' filter
    assert list(table.read) == ['r1', 'r3']
    sizes = table.set_index('read').repeat_size
    assert sizes['r1'] == pytest.approx(2000, rel=0.01)
    assert sizes['r3'] == pytest.approx(3000, rel=0.01)

    saved = pandas.read_csv(out_file, sep="\t")
    assert list(saved.read) == ['r1', 'r3']
    assert "running fft" in capsys.readouterr().out


def test_build_all_v_all_cmer_table_filters(sklearn, all_v_all):
    build = lambda **kw: concatemers.build_all_v_all_cmer_table(
        all_v_all, method=concatemers.CLUST, **kw)
    assert list(build(final_filter='all').read) == ['r1', 'r2', 'r3']
    assert list(build(final_filter='state').read) == ['r1', 'r3']
    # an unknown (non numeric) filter is the same as 'all'
    assert list(build(final_filter='whatever').read) == ['r1', 'r2', 'r3']


@pytest.mark.xfail(strict=True, raises=KeyError, reason=(
    "the numeric final_filter looks for a 'score' column, but the table "
    "column is 'repeat_score'"))
def test_build_all_v_all_cmer_table_score_filter(sklearn, all_v_all):
    table = concatemers.build_all_v_all_cmer_table(
        all_v_all, method=concatemers.CLUST, final_filter=0.5)
    assert list(table.read) == ['r1', 'r3']


def test_build_all_v_all_cmer_table_no_self_filter(sklearn, all_v_all):
    table = concatemers.build_all_v_all_cmer_table(
        all_v_all, method=concatemers.CLUST, final_filter='all',
        self_filter=False)
    # the cross hit is grouped with r1's hits (by query) but doesn't change
    # the result for this example
    assert list(table.read) == ['r1', 'r2', 'r3']


def test_build_all_v_all_cmer_table_from_file(scipy, all_v_all, tmp_path):
    path = tmp_path / "all_v_all.m8"
    path.write_text(to_m8(all_v_all))
    table = concatemers.build_all_v_all_cmer_table(str(path))
    assert list(table.read) == ['r1', 'r3']


def test_get_cmer_check_line(scipy, repeat_hits, tmp_path):
    out = tmp_path / "check.tsv"
    concatemers.get_cmer_check_line(repeat_hits, str(out), candidate='r1')
    fields = out.read_text().rstrip("\n").split("\t")
    assert len(fields) == 6
    assert fields[0] == 'r1'
    assert float(fields[1]) == pytest.approx(2000, rel=0.01)
    assert fields[4] == 'True'
    assert fields[5] == '10000'


@pytest.mark.parametrize("arg, value", [
    ("5", 5), ("-3", -3), ("2.5", 2.5), ("1e3", 1000.0),
    ("True", True), ("False", False), ("text", "text"), ("", ""),
])
def test_try_to_convert_arg(arg, value):
    converted = concatemers._try_to_convert_arg(arg)
    assert converted == value
    assert type(converted) == type(value)


def test_main(scipy, repeat_hits, tmp_path, monkeypatch, capsys):
    in_file = tmp_path / "self.m8"
    in_file.write_text(to_m8(repeat_hits))
    out_file = tmp_path / "check.tsv"
    monkeypatch.setattr(sys, 'argv', [
        'concatemers.py', 'get_cmer_check_line', str(in_file), str(out_file),
        'candidate=r1', 'qh_cutoff=500'])
    concatemers.main()
    fields = out_file.read_text().split("\t")
    assert fields[0] == 'r1'
    assert float(fields[1]) == pytest.approx(2000, rel=0.01)
    assert "running get_cmer_check_line" in capsys.readouterr().out


@pytest.mark.parametrize("argv, error, match", [
    ([], Exception, "function name"),
    (["os.system", "ls"], Exception, "word characters"),
    (["no_such_function"], NameError, "no_such_function"),
])
def test_main_errors(monkeypatch, argv, error, match):
    monkeypatch.setattr(sys, 'argv', ['concatemers.py'] + argv)
    with pytest.raises(error, match=match):
        concatemers.main()


###
# plotting


def test_dotplot(ax, repeat_hits):
    concatemers.dotplot(repeat_hits, ax=ax)
    lines = ax.get_lines()
    assert len(lines) == 9
    # scaled to the longest hit
    assert max(max(l.get_xdata()) for l in lines) == 1.0

    concatemers.dotplot(repeat_hits, ax=ax, seq_length=1)
    assert max(max(l.get_xdata()) for l in ax.get_lines()) == 10000


def test_dotplot_makes_axes(close_figures, repeat_hits):
    from matplotlib import pyplot as plt
    concatemers.dotplot(repeat_hits, subplots_kws={'figsize': [3, 3]})
    assert list(plt.gcf().get_size_inches()) == [3, 3]
