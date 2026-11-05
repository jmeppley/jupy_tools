import io

import numpy
import pandas
import pytest

from jme.jupy_tools import hit_tables


def _stream(text):
    return io.StringIO(text)


def _write(tmp_path, text, name="hits.tsv"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


###
# Test data

M8_TEXT = 'tp1:k4.min_d0.1.nn15:d004c256-f168-4f0c-9e42-03dacdc5fd30\ttp1:k4.min_d0.5.nn15:1a8c16f8-acb3-452e-83c1-4eeab66f5a85\t96.32\t11162\t358\t28\t11150\t1\t23325\t34445\t0\t1.6e+04\t34135\t34446\t10144\ntp1:k4.min_d0.1.nn15:d004c256-f168-4f0c-9e42-03dacdc5fd30\ttp1:k5.min_d0.5.nn15:1a8c16f8-acb3-452e-83c1-4eeab66f5a85\t96.31\t11162\t358\t29\t11150\t1\t23369\t34488\t0\t1.6e+04\t34135\t34491\t10135\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp1:k6.min_d0.nn20:335cb453-9678-48d9-a9df-8da26bc90aef\t98.73\t20901\t171\t64\t34928\t14078\t1\t20857\t0\t3.15e+04\t34944\t35130\t19923\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp1:k4.min_d0.5.nn10:ec5aceaf-f0c7-41cb-bdbd-82eb386082db\t98.72\t20901\t174\t64\t34928\t14078\t1\t20857\t0\t3.15e+04\t34944\t35118\t19917\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp1:k6.min_d0.5.nn5:ec5aceaf-f0c7-41cb-bdbd-82eb386082db\t98.68\t20900\t181\t63\t34928\t14078\t1\t20855\t0\t3.15e+04\t34944\t35096\t19909\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp2:k6.min_d0.nn5:4d6b63f3-c2f3-4d55-a5b1-06194b58ed44\t98.27\t20933\t236\t72\t34928\t14078\t1\t20889\t0\t3.12e+04\t34944\t35202\t19705\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp2:k6.min_d0.nn10:4d6b63f3-c2f3-4d55-a5b1-06194b58ed44\t98.29\t20936\t228\t75\t34928\t14078\t1\t20892\t0\t3.11e+04\t34944\t35196\t19697\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp2:k5.min_d0.1.nn5:4d6b63f3-c2f3-4d55-a5b1-06194b58ed44\t98.14\t20934\t263\t69\t34932\t14078\t1\t20887\t0\t3.11e+04\t34944\t35190\t19673\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp1:k6.min_d0.nn20:335cb453-9678-48d9-a9df-8da26bc90aef\t96.31\t10323\t326\t40\t11203\t901\t24000\t34287\t0\t1.47e+04\t34944\t35130\t9281\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp2:k6.min_d0.nn10:4d6b63f3-c2f3-4d55-a5b1-06194b58ed44\t96.56\t7961\t241\t29\t11203\t3259\t24032\t31975\t0\t1.14e+04\t34944\t35196\t7210\n'

PAF_TEXT = '437a607c-572a-4151-aced-2ecbaaaf2dd6\t80699\t781\t53821\t-\tcontig_18\t127941\t70093\t124878\t18878\t55497\t60\ttp:A:P\tcm:i:2404\ts1:i:18410\ts2:i:12981\tdv:f:0.0946\trl:i:747\n437a607c-572a-4151-aced-2ecbaaaf2dd6\t80699\t65401\t80609\t+\tcontig_326\t41785\t19929\t35247\t6566\t15929\t52\ttp:A:P\tcm:i:835\ts1:i:6373\ts2:i:5385\tdv:f:0.0792\trl:i:747\ne71fd889-4d97-4e3e-adad-52dbf7427f58\t75237\t55751\t75143\t+\tcontig_152\t117520\t12565\t32075\t8608\t19859\t19\ttp:A:P\tcm:i:1146\ts1:i:8492\ts2:i:8008\tdv:f:0.0760\trl:i:2049\ne71fd889-4d97-4e3e-adad-52dbf7427f58\t75237\t55679\t73934\t+\tcontig_1150\t56894\t9894\t28372\t8157\t18893\t0\ttp:A:S\tcm:i:1081\ts1:i:8008\tdv:f:0.0757\trl:i:2049\n7c41d647-ebf7-4dde-80e6-cfd18855ea93\t72268\t32537\t72194\t-\tcontig_2\t215798\t102920\t143071\t28961\t40578\t60\ttp:A:P\tcm:i:4176\ts1:i:28861\ts2:i:3414\tdv:f:0.0383\trl:i:275\n7c41d647-ebf7-4dde-80e6-cfd18855ea93\t72268\t14\t32485\t-\tcontig_577\t54763\t22091\t54755\t22826\t33066\t60\ttp:A:P\tcm:i:3302\ts1:i:22757\ts2:i:369\tdv:f:0.0409\trl:i:275\nd5b097c7-3418-4a47-b6d1-89abe4685b92\t70279\t50616\t69419\t+\tcontig_2284\t22779\t16\t19816\t7162\t19921\t0\ttp:A:P\tcm:i:932\ts1:i:6981\ts2:i:6968\tdv:f:0.0880\trl:i:201\nd5b097c7-3418-4a47-b6d1-89abe4685b92\t70279\t47390\t67102\t+\tcontig_854\t61810\t35692\t55802\t7057\t20250\t0\ttp:A:S\tcm:i:902\ts1:i:6968\tdv:f:0.0934\trl:i:201\n7e99fca7-47f8-4f8c-b060-ab3f38115fd4\t65629\t11877\t35094\t-\tcontig_1013\t89882\t36694\t60618\t2601\t24321\t60\ttp:A:P\tcm:i:289\ts1:i:2357\ts2:i:1056\tdv:f:0.1792\trl:i:2853\n7e99fca7-47f8-4f8c-b060-ab3f38115fd4\t65629\t42349\t65584\t-\tcontig_875\t161459\t100414\t124088\t2495\t24100\t60\ttp:A:P\tcm:i:281\ts1:i:2277\ts2:i:1300\tdv:f:0.1812\trl:i:2853\n5eb578b0-94df-436d-8b15-1c530394214c\t64108\t4689\t53038\t+\tcontig_1253\t108969\t59766\t108945\t23193\t49677\t60\ttp:A:P\tcm:i:3053\ts1:i:22990\ts2:i:14250\tdv:f:0.0720\trl:i:430\n5eb578b0-94df-436d-8b15-1c530394214c\t64108\t34656\t63995\t+\tcontig_341\t38096\t6187\t35393\t6081\t30031\t60\ttp:A:P\tcm:i:733\ts1:i:5840\ts2:i:1761\tdv:f:0.1332\trl:i:430\nb57d5447-df17-4619-abcd-48ea454f1eea\t63531\t3715\t52749\t+\tcontig_468\t121268\t63518\t114205\t14813\t51313\t60\ttp:A:P\tcm:i:1799\ts1:i:14354\ts2:i:6413\tdv:f:0.1063\trl:i:3871\nb57d5447-df17-4619-abcd-48ea454f1eea\t63531\t36915\t63424\t+\tcontig_557\t28416\t1619\t28327\t2151\t27315\t60\ttp:A:P\tcm:i:232\ts1:i:1908\ts2:i:1191\tdv:f:0.2010\trl:i:3871'

# A small, hand-checkable BLAST+ table (with qlen, hlen, raw_score columns)
#   q1/h1 has 3 fragments: the third overlaps the first on the query
#   q1/h2 is on the reverse strand of the hit
SMALL_BLAST_PLUS = (
    "# a comment line\n"
    #query hit  pctid mlen mm gaps qstart qend hstart hend evalue score qlen hlen raw
    "q1\th1\t90.0\t100\t10\t0\t1\t100\t1\t100\t1e-50\t200\t1000\t2000\t150\n"
    "q1\th1\t80.0\t50\t10\t0\t201\t250\t301\t350\t1e-20\t90\t1000\t2000\t70\n"
    "q1\th1\t70.0\t60\t18\t0\t50\t109\t500\t559\t1e-10\t60\t1000\t2000\t50\n"
    "q1\th2\t95.0\t200\t10\t0\t300\t499\t1000\t801\t1e-80\t350\t1000\t1500\t300\n"
    "q2\th1\t99.0\t300\t3\t0\t1\t300\t1\t300\t1e-90\t500\t400\t2000\t450\n"
)
# same thing without the 3 extra columns
SMALL_BLAST = "".join(
    "\t".join(line.split("\t")[:12]) + "\n"
    for line in SMALL_BLAST_PLUS.splitlines()
    if not line.startswith("#")
)

# lastal tabular output (hit, aka reference, comes before query)
LASTAL_TEXT = (
    "# LAST version 1234\n"
    #score hit start alnsize strand seqsize query start alnsize strand seqsize blocks EG2 E
    "150\th1\t10\t100\t+\t2000\tq1\t0\t102\t+\t1000\t50,2:0,48\tEG2=1e-30\tE=2.5e-40\n"
    "80\th2\t500\t50\t+\t1500\tq1\t400\t50\t-\t1000\t20,0:3,27\tEG2=1e-10\tE=3e-15\n"
)


###
# constants


def test_column_definitions():
    assert hit_tables.HIT_TABLE_COLUMNS[hit_tables.BLAST] == hit_tables.BLAST_COLUMNS
    assert len(hit_tables.BLAST_COLUMNS) == 12
    assert hit_tables.BLAST_PLUS_COLUMNS[:12] == hit_tables.BLAST_COLUMNS
    assert len(hit_tables.BLAST_PLUS_COLUMNS) == 15
    assert len(hit_tables.PAF_COLUMNS) == 12
    assert len(hit_tables.HIT_TABLE_COLUMNS[hit_tables.PAF_ALL]) == 16
    assert len(hit_tables.LASTAL_COLUMNS) == 14


@pytest.mark.parametrize("fmt, text", [
    (hit_tables.BLAST, SMALL_BLAST),
    (hit_tables.BLAST_PLUS, SMALL_BLAST_PLUS.splitlines()[1] + "\n"),
    (hit_tables.PAF, PAF_TEXT.splitlines()[0] + "\n"),
])
def test_hit_table_rexps(fmt, text):
    m = hit_tables.HIT_TABLE_REXPS[fmt].search(text)
    assert m is not None
    data = m.groupdict()
    for key in ['query', 'hit', 'qstart', 'qend', 'hstart', 'hend']:
        assert key in data


###
# computeLastHitValues


@pytest.mark.parametrize("blocks, mlen", [
    ("100", 100),
    ("73,0:1,15,0:1,13,0:1,9", 73 + 1 + 15 + 1 + 13 + 1 + 9),
    ("10,5:0,10", 25),
    ("10,2:7,10", 27),
])
def test_compute_last_hit_values(blocks, mlen):
    assert hit_tables.computeLastHitValues(blocks) == mlen


###
# parse_blast_m8 (pandas)


def test_parse_m8():
    df = hit_tables.parse_blast_m8(_stream(M8_TEXT), format=hit_tables.BLAST_PLUS)
    assert df.shape == (10, 15)
    assert len(set(df.hit)) == 8
    assert df.iloc[0, 12] == 34135
    assert df.sort_values('pctid').pctid[0] == 96.32


def test_parse_paf():
    df = hit_tables.parse_blast_m8(_stream(PAF_TEXT), format=hit_tables.PAF)
    assert df.shape == (14, 13)
    assert len(set(df.hit)) == 14
    assert len(set(df['query'])) == 7
    assert df.iloc[0, 10] == 55497
    assert numpy.round(df.sort_values('pctid').pctid[0], 3) == 34.016


def test_parse_paf_pctid():
    df = hit_tables.parse_blast_m8(_stream(PAF_TEXT), format=hit_tables.PAF)
    assert numpy.allclose(df.pctid, 100 * df.matches / df.mlen)
    assert list(df.columns) == hit_tables.PAF_COLUMNS + ['pctid']


def test_parse_paf_all():
    df = hit_tables.parse_blast_m8(_stream(PAF_TEXT), format=hit_tables.PAF_ALL)
    assert df.shape == (14, 17)
    assert df.tp.iloc[0] == 'tp:A:P'
    assert df.cm.iloc[0] == 'cm:i:2404'


@pytest.mark.xfail(strict=True, reason=(
    "PAF_ALL names tag columns by position (tp, cm, dv, rl), but minimap2 "
    "writes s1 (and usually s2) between cm and dv"))
def test_parse_paf_all_tag_columns():
    df = hit_tables.parse_blast_m8(_stream(PAF_TEXT), format=hit_tables.PAF_ALL)
    assert df.dv.iloc[0] == 'dv:f:0.0946'
    assert df.rl.iloc[0] == 'rl:i:747'


def test_parse_blast():
    df = hit_tables.parse_blast_m8(_stream(SMALL_BLAST))
    assert list(df.columns) == hit_tables.BLAST_COLUMNS
    assert df.shape == (5, 12)
    assert df.evalue.iloc[0] == 1e-50


def test_parse_blast_from_path(tmp_path):
    path = _write(tmp_path, SMALL_BLAST_PLUS)
    df = hit_tables.parse_blast_m8(path, format=hit_tables.BLAST_PLUS)
    # comment line is skipped
    assert df.shape == (5, 15)
    assert list(df['query']) == ['q1', 'q1', 'q1', 'q1', 'q2']


def test_parse_blast_extra_columns_ignored():
    """ BLAST format only uses the first 12 columns """
    df = hit_tables.parse_blast_m8(_stream(SMALL_BLAST_PLUS))
    assert df.shape == (5, 12)


def test_parse_blast_skiprows():
    text = "garbage header\n" + SMALL_BLAST
    df = hit_tables.parse_blast_m8(_stream(text), skiprows=1)
    assert df.shape == (5, 12)


def test_parse_blast_cutoffs():
    parse = lambda **c: hit_tables.parse_blast_m8(
        _stream(SMALL_BLAST_PLUS), format=hit_tables.BLAST_PLUS, **c)

    # positive cutoff: minimum value
    assert sorted(parse(pctid=90).pctid) == [90, 95, 99]
    # combined
    assert list(parse(pctid=90, mlen=150).mlen) == [200, 300]
    # negative cutoff: maximum value
    assert sorted(parse(evalue=-1e-30).evalue) == [1e-90, 1e-80, 1e-50]
    # zero is ignored
    assert parse(pctid=0).shape[0] == 5


def test_parse_blast_empty_returns_none():
    assert hit_tables.parse_blast_m8(_stream("# nothing but comments\n")) is None


def test_parse_blast_cutoffs_remove_all():
    df = hit_tables.parse_blast_m8(_stream(SMALL_BLAST), pctid=100)
    assert df.shape[0] == 0


def test_parse_lastal():
    df = hit_tables.parse_blast_m8(_stream(LASTAL_TEXT), format=hit_tables.LASTAL)
    assert df.shape[0] == 2
    assert list(df.hit) == ['h1', 'h2']
    assert list(df['query']) == ['q1', 'q1']
    # match lengths computed from blocks
    assert list(df.mlen) == [100, 50]
    # E= prefix removed
    assert list(df.evalue) == [2.5e-40, 3e-15]
    # end positions computed from start + alignment length
    assert list(df.hend) == [109, 549]
    assert list(df.qend) == [101, 351]
    assert list(df.qsmult) == [1, -1]


def test_parse_blast_pair_overlap_buffer():
    df = hit_tables.parse_blast_m8(_stream(SMALL_BLAST_PLUS),
                                   format=hit_tables.BLAST_PLUS,
                                   pair_overlap_buffer=0)
    # the 3rd fragment of q1/h1 overlaps the 1st on the query
    assert df.shape[0] == 4
    assert 70.0 not in set(df.pctid)

    # with cutoffs applied after overlap removal
    df = hit_tables.parse_blast_m8(_stream(SMALL_BLAST_PLUS),
                                   format=hit_tables.BLAST_PLUS,
                                   pair_overlap_buffer=0, pctid=85)
    assert sorted(df.pctid) == [90, 95, 99]


def test_parse_blast_pair_overlap_bad_format():
    with pytest.raises(Exception, match="not supported"):
        hit_tables.parse_blast_m8(_stream(LASTAL_TEXT),
                                  format=hit_tables.LASTAL,
                                  pair_overlap_buffer=0)


###
# parse_blast_m8 (polars)


@pytest.fixture
def pl():
    return pytest.importorskip("polars")


def test_parse_blast_polars(pl, tmp_path):
    path = _write(tmp_path, SMALL_BLAST_PLUS)
    df = hit_tables.parse_blast_m8(path, format=hit_tables.BLAST_PLUS,
                                   use_polars=True)
    assert isinstance(df, pl.DataFrame)
    assert df.shape == (5, 15)
    assert df.columns == hit_tables.BLAST_PLUS_COLUMNS
    assert df['qlen'].to_list() == [1000, 1000, 1000, 1000, 400]


def test_parse_blast_polars_lazy(pl, tmp_path):
    path = _write(tmp_path, SMALL_BLAST)
    lf = hit_tables.parse_blast_m8(path, use_polars="lazy")
    assert isinstance(lf, pl.LazyFrame)
    assert lf.collect().shape == (5, 12)

    # any other string is not lazy
    df = hit_tables.parse_blast_m8(path, use_polars="eager")
    assert isinstance(df, pl.DataFrame)


def test_parse_blast_polars_cutoffs(pl, tmp_path):
    path = _write(tmp_path, SMALL_BLAST_PLUS)
    df = hit_tables.parse_blast_m8(path, format=hit_tables.BLAST_PLUS,
                                   use_polars=True, pctid=90, evalue=-1e-60,
                                   mlen=0)
    assert sorted(df['pctid'].to_list()) == [95, 99]


def test_parse_blast_polars_empty(pl, tmp_path):
    path = _write(tmp_path, SMALL_BLAST)
    assert hit_tables.parse_blast_m8(path, use_polars=True, pctid=100) is None


def test_parse_paf_polars(pl, tmp_path):
    # just the 12 standard columns (see the xfail below for tags)
    path = _write(tmp_path, "".join(
        "\t".join(line.split("\t")[:12]) + "\n"
        for line in PAF_TEXT.splitlines()
    ))
    df = hit_tables.parse_blast_m8(path, format=hit_tables.PAF, use_polars=True)
    pd_df = hit_tables.parse_blast_m8(path, format=hit_tables.PAF)
    assert df.shape == pd_df.shape
    assert numpy.allclose(df['pctid'].to_numpy(), pd_df.pctid.values)


@pytest.mark.xfail(strict=True, reason=(
    "the polars parser passes new_columns to scan_csv, which fails when the "
    "file has more columns than the format (eg minimap2 PAF tags)"))
def test_parse_paf_polars_with_tags(pl, tmp_path):
    path = _write(tmp_path, PAF_TEXT)
    df = hit_tables.parse_blast_m8(path, format=hit_tables.PAF, use_polars=True)
    assert df.shape == (14, 13)


def test_parse_blast_polars_matches_pandas(pl, tmp_path):
    path = _write(tmp_path, M8_TEXT)
    df = hit_tables.parse_blast_m8(path, format=hit_tables.BLAST_PLUS,
                                   use_polars=True)
    pd_df = hit_tables.parse_blast_m8(path, format=hit_tables.BLAST_PLUS)
    assert df.shape == pd_df.shape
    assert df['hit'].to_list() == list(pd_df.hit)
    assert df['score'].to_list() == list(pd_df.score)


@pytest.mark.xfail(strict=True, raises=NameError, reason=(
    "_parse_blast_m8_polars calls return_dtype(...) instead of passing "
    "return_dtype=... for the LASTAL qstrand column"))
def test_parse_lastal_polars(pl, tmp_path):
    path = _write(tmp_path, LASTAL_TEXT)
    df = hit_tables.parse_blast_m8(path, format=hit_tables.LASTAL,
                                   use_polars=True)
    pd_df = hit_tables.parse_blast_m8(path, format=hit_tables.LASTAL)
    assert df['mlen'].to_list() == list(pd_df.mlen)
    assert df['hend'].to_list() == list(pd_df.hend)
    assert df['qend'].to_list() == list(pd_df.qend)


###
# overlap removal


def test_generate_nonoverlapping_lines():
    lines = list(hit_tables.generate_nonoverlapping_lines(
        _stream(SMALL_BLAST_PLUS), format=hit_tables.BLAST_PLUS))
    assert len(lines) == 4
    assert not any(l.startswith("#") for l in lines)
    assert all(l.endswith("\n") for l in lines)
    assert not any("\t70.0\t" in l for l in lines)


def test_generate_nonoverlapping_lines_from_path(tmp_path):
    path = _write(tmp_path, SMALL_BLAST)
    lines = list(hit_tables.generate_nonoverlapping_lines(
        path, format=hit_tables.BLAST))
    assert len(lines) == 4


def test_generate_nonoverlapping_lines_buffer():
    text = (
        # second hit overlaps first by 10bp on both query and hit
        "q\th\t90\t100\t0\t0\t1\t100\t1\t100\t0\t100\n"
        "q\th\t90\t100\t0\t0\t91\t190\t91\t190\t0\t100\n"
        # same coordinates, but a different pair: always kept
        "q\tother\t90\t100\t0\t0\t1\t100\t1\t100\t0\t100\n"
        # a short hit (shorter than the buffer) inside the first: kept
        # only if the buffer is bigger than its length
        "q\th\t90\t5\t0\t0\t10\t14\t10\t14\t0\t10\n"
    )
    keep = lambda b: [
        l.split("\t")[:2] + l.split("\t")[6:8]
        for l in hit_tables.generate_nonoverlapping_lines(
            _stream(text), format=hit_tables.BLAST, buffer=b)
    ]
    assert keep(0) == [['q', 'h', '1', '100'], ['q', 'other', '1', '100']]
    assert keep(10) == [['q', 'h', '1', '100'], ['q', 'h', '91', '190'],
                        ['q', 'other', '1', '100'], ['q', 'h', '10', '14']]


def test_generate_nonoverlapping_lines_reverse_strand():
    text = (
        "q\th\t90\t100\t0\t0\t1\t100\t200\t101\t0\t100\n"
        # overlaps the first hit on the (reversed) hit coordinates
        "q\th\t90\t100\t0\t0\t301\t400\t150\t110\t0\t100\n"
        "q\th\t90\t100\t0\t0\t501\t600\t300\t201\t0\t100\n"
    )
    lines = list(hit_tables.generate_nonoverlapping_lines(
        _stream(text), format=hit_tables.BLAST))
    assert [l.split("\t")[6] for l in lines] == ['1', '501']


def test_generate_nonoverlapping_lines_paf():
    # the same alignment twice
    text = PAF_TEXT.splitlines()[0] + "\n"
    lines = list(hit_tables.generate_nonoverlapping_lines(
        _stream(text * 2), format=hit_tables.PAF))
    assert len(lines) == 1


def test_generate_nonoverlapping_lines_bad_line(capsys):
    with pytest.raises(AttributeError):
        list(hit_tables.generate_nonoverlapping_lines(
            _stream("not\ta\thit\tline\n"), format=hit_tables.BLAST))
    assert "ERROR on line" in capsys.readouterr().out


def test_remove_pair_overlaps():
    with hit_tables.remove_pair_overlaps(_stream(SMALL_BLAST),
                                         format=hit_tables.BLAST) as handle:
        content = handle.read().decode()
    assert len(content.splitlines()) == 4


###
# aggregation (pandas)


def _check_small_agg(agg, keep_extent=False):
    """ shared checks for aggregated SMALL_BLAST_PLUS (w/ overlap removal) """
    agg = agg.set_index(['query', 'hit'])
    assert len(agg) == 3

    q1h1 = agg.loc[('q1', 'h1')]
    assert q1h1.matches == pytest.approx(130)
    assert q1h1.mlen == 150
    assert q1h1.pctid == pytest.approx(100 * 130 / 150)
    assert q1h1.qmlen == 150
    assert q1h1.hmlen == 150
    assert q1h1.qlen == 1000
    assert q1h1.hlen == 2000
    assert q1h1.mfrac == pytest.approx(200 * 130 / 3000)

    q1h2 = agg.loc[('q1', 'h2')]
    assert q1h2.hmlen == 200
    assert q1h2.pctid == pytest.approx(95)

    if keep_extent:
        assert (q1h1.qstart, q1h1.qend, q1h1.hstart, q1h1.hend) == \
            (1, 250, 1, 350)
        # reverse strand coordinates are sorted
        assert (q1h2.hstart, q1h2.hend) == (801, 1000)


def test_agg_hit_df():
    hits = hit_tables.parse_blast_m8(_stream(SMALL_BLAST_PLUS),
                                     format=hit_tables.BLAST_PLUS,
                                     pair_overlap_buffer=0)
    agg = hit_tables.agg_hit_df(hits)
    _check_small_agg(agg)
    assert 'qstart' not in agg.columns


def test_agg_hit_df_keep_extent():
    hits = hit_tables.parse_blast_m8(_stream(SMALL_BLAST_PLUS),
                                     format=hit_tables.BLAST_PLUS,
                                     pair_overlap_buffer=0)
    agg = hit_tables.agg_hit_df(hits, keep_extent=True)
    _check_small_agg(agg, keep_extent=True)


def test_agg_hit_df_max_fragments():
    hits = hit_tables.parse_blast_m8(_stream(SMALL_BLAST_PLUS),
                                     format=hit_tables.BLAST_PLUS,
                                     pair_overlap_buffer=0)
    agg = hit_tables.agg_hit_df(hits, max_hit_fragments=1) \
        .set_index(['query', 'hit'])
    # only the best fragment of q1/h1 is counted
    assert agg.loc[('q1', 'h1'), 'matches'] == pytest.approx(90)
    assert agg.loc[('q1', 'h1'), 'mlen'] == 100
    assert agg.loc[('q1', 'h1'), 'pctid'] == pytest.approx(90)


def test_agg_hit_df_without_lengths():
    """ plain BLAST tables have no qlen/hlen, so no mfrac """
    hits = hit_tables.parse_blast_m8(_stream(SMALL_BLAST))
    agg = hit_tables.agg_hit_df(hits)
    assert 'mfrac' not in agg.columns
    assert 'qlen' not in agg.columns
    # no overlap removal, so all 3 q1/h1 fragments are summed
    q1h1 = agg.set_index(['query', 'hit']).loc[('q1', 'h1')]
    assert q1h1.mlen == 210
    assert q1h1.matches == pytest.approx(90 + 40 + 42)


def test_agg_hit_df_with_matches_column():
    """ PAF tables have a matches column, so it's used directly """
    hits = hit_tables.parse_blast_m8(_stream(PAF_TEXT), format=hit_tables.PAF)
    agg = hit_tables.agg_hit_df(hits)
    assert len(agg) == 14
    assert agg.matches.sum() == hits.matches.sum()
    assert 'mfrac' in agg.columns


def test_agg_hit_table():
    agg = hit_tables.agg_hit_table(_stream(SMALL_BLAST_PLUS),
                                   format=hit_tables.BLAST_PLUS)
    _check_small_agg(agg)


def test_agg_hit_table_ovl_buffer():
    """ a big enough buffer lets the overlapping fragment through """
    agg = hit_tables.agg_hit_table(_stream(SMALL_BLAST_PLUS),
                                   format=hit_tables.BLAST_PLUS,
                                   ovl_buffer=100)
    q1h1 = agg.set_index(['query', 'hit']).loc[('q1', 'h1')]
    assert q1h1.mlen == 210


def test_agg_hit_table_empty():
    assert hit_tables.agg_hit_table(_stream("# empty\n")) is None


###
# aggregation (polars)


def _polars_to_pandas(df):
    return pandas.DataFrame(df.to_dict(as_series=False))


def test_agg_hit_table_polars(pl, tmp_path):
    path = _write(tmp_path, SMALL_BLAST_PLUS)
    agg = hit_tables.agg_hit_table(path, format=hit_tables.BLAST_PLUS,
                                   use_polars=True)
    assert isinstance(agg, pl.DataFrame)
    _check_small_agg(_polars_to_pandas(agg))


def test_agg_hit_table_polars_lazy(pl, tmp_path):
    path = _write(tmp_path, SMALL_BLAST_PLUS)
    agg = hit_tables.agg_hit_table(path, format=hit_tables.BLAST_PLUS,
                                   use_polars='lazy')
    assert isinstance(agg, pl.LazyFrame)
    _check_small_agg(_polars_to_pandas(agg.collect()))


@pytest.mark.xfail(strict=True, reason=(
    "agg_hit_df_polars names both the min and max extent expressions after "
    "the start column, so keep_extent=True raises a DuplicateError"))
def test_agg_hit_table_polars_keep_extent(pl, tmp_path):
    path = _write(tmp_path, SMALL_BLAST_PLUS)
    agg = hit_tables.agg_hit_table(path, format=hit_tables.BLAST_PLUS,
                                   use_polars=True, keep_extent=True)
    _check_small_agg(_polars_to_pandas(agg), keep_extent=True)


def test_agg_hit_df_polars_max_fragments(pl, tmp_path):
    path = _write(tmp_path, SMALL_BLAST_PLUS)
    hits = hit_tables.parse_blast_m8(path, format=hit_tables.BLAST_PLUS,
                                     use_polars=True, pair_overlap_buffer=0)
    agg = _polars_to_pandas(
        hit_tables.agg_hit_df_polars(hits, max_hit_fragments=1)
    ).set_index(['query', 'hit'])
    assert agg.loc[('q1', 'h1'), 'matches'] == pytest.approx(90)
    assert agg.loc[('q1', 'h1'), 'mlen'] == 100


def test_agg_hit_df_polars_without_lengths(pl, tmp_path):
    path = _write(tmp_path, SMALL_BLAST)
    hits = hit_tables.parse_blast_m8(path, use_polars=True)
    agg = hit_tables.agg_hit_df_polars(hits)
    assert 'mfrac' not in agg.columns
    assert 'pctid' in agg.columns


def test_agg_hit_df_dispatches_to_polars(pl, tmp_path):
    """ agg_hit_df hands non-pandas tables to agg_hit_polars_df """
    path = _write(tmp_path, SMALL_BLAST_PLUS)
    hits = hit_tables.parse_blast_m8(path, format=hit_tables.BLAST_PLUS,
                                     use_polars=True, pair_overlap_buffer=0)
    agg = hit_tables.agg_hit_df(hits)
    assert isinstance(agg, pl.DataFrame)
    _check_small_agg(_polars_to_pandas(agg))

    lazy_agg = hit_tables.agg_hit_df(hits.lazy(), keep_extent=True)
    assert isinstance(lazy_agg, pl.LazyFrame)
    _check_small_agg(_polars_to_pandas(lazy_agg.collect()), keep_extent=True)


def test_agg_hit_polars_df_errors(pl):
    with pytest.raises(Exception, match="can't handle type"):
        hit_tables.agg_hit_polars_df([1, 2, 3])

    hits = pl.DataFrame({'query': ['q'], 'hit': ['h'], 'pctid': [90.0],
                         'mlen': [100]})
    with pytest.raises(Exception, match="max_hit_fragments"):
        hit_tables.agg_hit_polars_df(hits, max_hit_fragments=2)


def test_check_columns(pl):
    df = pl.DataFrame({'a': [1], 'b': [2]})
    assert hit_tables._check_columns(df, 'a')
    assert not hit_tables._check_columns(df, 'c')
    assert hit_tables._check_columns(df.lazy(), 'b', is_lazy=True)
    assert not hit_tables._check_columns(df.lazy(), 'c', is_lazy=True)


###
# SAM flags


def test_decode_sam_flags_zero():
    flags = hit_tables.decode_sam_flags(0)
    assert set(flags) == set(hit_tables.bit_names)
    assert not any(flags.values())


@pytest.mark.parametrize("flag, expected", [
    (1, {'paired'}),
    (4, {'unmapped'}),
    (16, {'read rev strand'}),
    (99, {'paired', 'mapped in pair', 'mate rev strand', 'first in pair'}),
    (147, {'paired', 'mapped in pair', 'read rev strand', 'second in pair'}),
    (256, {'not primary alig'}),
    (2048, {'supplementary'}),
    (4095, set(hit_tables.bit_names)),
    ("83", {'paired', 'mapped in pair', 'read rev strand', 'first in pair'}),
])
def test_decode_sam_flags(flag, expected):
    flags = hit_tables.decode_sam_flags(flag)
    assert {k for k, v in flags.items() if v} == expected


@pytest.mark.parametrize("flag, primary", [
    (0, True), (16, True), (99, True), (147, True),
    (4, False), (256, False), (2048, False), (2064, False),
])
def test_is_sam_flag_primary_hit(flag, primary):
    assert hit_tables.is_sam_flag_primary_hit(flag) == primary


@pytest.mark.parametrize("flag, first", [
    (99, True), (83, True), (147, False), (0, False),
])
def test_is_sam_flag_first_in_pair(flag, first):
    assert hit_tables.is_sam_flag_first_in_pair(flag) == first
