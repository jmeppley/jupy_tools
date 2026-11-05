import io
import logging
import os

import numpy
import pandas
import pytest

from jme.jupy_tools import utils


###
# LogLogger


def test_log_logger(tmp_path):
    """ messages go to the logger with exponentially decreasing frequency """
    test_log_file = str(tmp_path / 'logger.test.out')
    logger = logging.getLogger('test_log_logger')
    file_handler = logging.FileHandler(test_log_file)
    logger.addHandler(file_handler)
    logger.propagate = False
    try:
        logger.warning("Init logger")
        llogger = utils.LogLogger(3, 1, 2, logger=logger, level=logging.WARNING)
        for i in range(10):
            llogger.log(f"Iteration: {i}")
    finally:
        logger.removeHandler(file_handler)
        file_handler.close()

    with open(test_log_file) as lines:
        lines = list(lines)
    assert len(lines) == 3
    assert lines[0] == 'Init logger\n'
    assert lines[1] == 'Iteration: 0\n'
    assert lines[2] == 'Iteration: 4\n'


def test_log_logger_print(capsys):
    """ with no logger, output is printed, and %-style args are formatted lazily """
    llogger = utils.LogLogger()
    for i in range(20):
        llogger.log("iteration: %s", i)
    printed = capsys.readouterr().out.splitlines()
    # matches the example in the class docstring
    assert printed == [f"iteration: {i}" for i in [0, 1, 3, 7, 15]]


def test_log_logger_no_format_args(capsys):
    llogger = utils.LogLogger()
    llogger.log("just a message")
    assert capsys.readouterr().out == "just a message\n"


def test_log_logger_check_pos():
    llogger = utils.LogLogger(multiplier=10)
    triggers = [i for i in range(1000) if llogger.check_pos()]
    assert triggers == [0, 9, 99, 999]


@pytest.mark.parametrize("level, method", [
    (logging.DEBUG, 'debug'),
    (logging.DEBUG - 1, 'debug'),
    (logging.INFO, 'info'),
    (logging.WARNING, 'warning'),
    (logging.ERROR, 'error'),
    (logging.CRITICAL, 'critical'),
])
def test_get_print_func(level, method):
    logger = logging.getLogger('test_get_print_func')
    assert utils._get_print_func(logger, level) == getattr(logger, method)


def test_get_print_func_no_logger():
    assert utils._get_print_func(None, logging.INFO) is utils._print


###
# small helpers


def test_first():
    assert utils.first([3, 2, 1]) == 3
    assert utils.first("abc") == "a"
    assert utils.first({'only': 1}) == 'only'
    assert utils.first(i * 2 for i in range(5, 10)) == 10
    with pytest.raises(StopIteration):
        utils.first([])


def test_classify():
    classes = utils.classify(range(10), lambda i: i % 3)
    assert dict(classes) == {0: [0, 3, 6, 9], 1: [1, 4, 7], 2: [2, 5, 8]}
    # it's a defaultdict, so missing classes are empty lists
    assert classes['missing'] == []


def test_classify_empty():
    assert dict(utils.classify([], len)) == {}


def test_read_tsv():
    tsv = io.StringIO("name\ta\tb\nx\t1\t2\ny\t3\t4\n")
    df = utils.read_tsv(tsv)
    assert list(df.index) == ['x', 'y']
    assert list(df.columns) == ['a', 'b']
    assert df.loc['y', 'b'] == 4


def test_read_tsv_kwargs_override_defaults():
    csv = io.StringIO("name,a,b\nx,1,2\n")
    df = utils.read_tsv(csv, sep=',', index_col=None)
    assert list(df.columns) == ['name', 'a', 'b']
    assert df.shape == (1, 3)


###
# best fit


def test_best_fit():
    x = numpy.array([0, 1, 2, 3, 4])
    for i in range(1, 10, 2):
        y = i * x
        slope, intercept = utils.get_best_fit(x, y, True)
        assert abs(slope - i) < 1e-10
        assert intercept == 0

        slope, intercept = utils.get_best_fit(x, y)
        assert abs(slope - i) < 1e-10
        assert abs(intercept) < 1e-10


def test_best_fit_with_intercept():
    x = numpy.array([0, 1, 2, 3, 4])
    y = 2 * x + 5
    slope, intercept = utils.get_best_fit(x, y)
    assert slope == pytest.approx(2)
    assert intercept == pytest.approx(5)

    # forcing the line through zero gives a steeper slope
    slope, intercept = utils.get_best_fit(x, y, force_intercept=True)
    assert intercept == 0
    assert slope > 2


def test_best_fit_force_slope():
    x = numpy.array([0, 1, 2, 3, 4])
    y = x + 3
    slope, intercept = utils.get_best_fit(x, y, force_slope=True)
    assert slope == 1
    assert intercept == pytest.approx(3)


def test_best_fit_force_both():
    x = numpy.array([0, 1, 2, 3, 4])
    assert tuple(utils.get_best_fit(x, 7 * x, True, True)) == (1, 0)


###
# streams and commands


def test_iterable_to_stream():
    lines = ["a\tb\n", "1\t2\n", "3\t4\n"]
    stream = utils.iterable_to_stream(iter(lines))
    assert stream.read() == b"".join(l.encode() for l in lines)


def test_iterable_to_stream_bytes():
    chunks = [b"abc", b"def"]
    stream = utils.iterable_to_stream(iter(chunks), str_to_bytes=False)
    assert stream.read() == b"abcdef"


def test_iterable_to_stream_small_buffer():
    """ chunks bigger than the buffer get split across reads """
    data = ["x" * 50, "y" * 50]
    stream = utils.iterable_to_stream(iter(data), buffer_size=8)
    assert stream.read() == ("x" * 50 + "y" * 50).encode()


def test_iterable_to_stream_into_pandas():
    lines = ("{}\t{}\n".format(i, i * i) for i in range(5))
    df = pandas.read_csv(utils.iterable_to_stream(lines), sep="\t",
                         header=None, names=['i', 'sq'])
    assert list(df.sq) == [0, 1, 4, 9, 16]


def test_get_dataframe_from_cmd():
    df = utils.get_dataframe_from_cmd("printf 'a\\tb\\n1\\t2\\n3\\t4\\n'")
    assert list(df.columns) == ['a', 'b']
    assert list(df.b) == [2, 4]


def test_get_dataframe_from_cmd_kwargs():
    df = utils.get_dataframe_from_cmd(['printf', 'x,y\n5,6\n'], shell=False,
                                      sep=',', index_col=0)
    assert df.loc[5, 'y'] == 6


###
# parsers


def test_parse_prinseq_stats():
    stats = io.StringIO(
        "stats_info\treads\t1000\n"
        "stats_info\tbases\t150000\n"
        "stats_len\tmean\t150.0\n"
    )
    series = utils.parse_prinseq_stats(stats)
    assert isinstance(series, pandas.Series)
    assert list(series.index) == ['reads', 'bases', 'mean']
    assert series['reads'] == 1000
    assert series['mean'] == 150.0


EGGNOG_HEADER = (
    "## Thu Jan  1 00:00:00 2021\n"
    "## emapper-2.1.2\n"
    "## /path/to/emapper.py -i x.faa\n"
    "##\n"
    "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tDescription\n"
)
EGGNOG_ROWS = (
    "gene1\t123.abc\t1e-50\t200.0\tCOG0001@1|root,COG0001@2|Bacteria,4NEW1@1224|Proteobacteria\tsome protein\n"
    "gene2\t456.def\t1e-20\t90.5\tCOG0002@1|root,arCOG0002@2157|Archaea\tother protein\n"
    "gene3\t789.ghi\t1e-10\t50.1\tKOG0003@2759|Eukaryota\tyet another\n"
)
EGGNOG_FOOTER = (
    "## 3 queries scanned\n"
    "## Total time (seconds): 1.0\n"
)
EGGNOG_TEXT = EGGNOG_HEADER + EGGNOG_ROWS + EGGNOG_FOOTER


@pytest.fixture
def eggnog_file(tmp_path):
    path = tmp_path / "test.emapper.annotations"
    path.write_text(EGGNOG_TEXT)
    return str(path)


def test_parse_eggnog_annotations(eggnog_file):
    df = utils.parse_eggnog_annotations(eggnog_file)
    assert df.shape[0] == 3
    # leading "#" is stripped from first column name, which becomes the index
    assert df.index.name == 'query'
    assert list(df.index) == ['gene1', 'gene2', 'gene3']
    assert df.loc['gene2', 'score'] == 90.5
    assert list(df.root_og) == ['COG0001', 'COG0002', 'KOG0003']
    assert list(df.tax_og) == ['4NEW1', 'arCOG0002', 'KOG0003']


def test_parse_eggnog_annotations_no_ogs(eggnog_file):
    df = utils.parse_eggnog_annotations(eggnog_file, ogs_col=None)
    assert df.shape == (3, 5)
    assert 'root_og' not in df.columns


def test_parse_eggnog_annotations_polars(eggnog_file):
    pl = pytest.importorskip("polars")
    lazy = utils.parse_eggnog_annotations_polars(eggnog_file)
    assert isinstance(lazy, pl.LazyFrame)

    df = utils.parse_eggnog_annotations_polars(eggnog_file, lazy=False)
    assert isinstance(df, pl.DataFrame)
    assert df.height == 3
    assert df['query'].to_list() == ['gene1', 'gene2', 'gene3']
    assert df['root_og'].to_list() == ['COG0001', 'COG0002', 'KOG0003']
    assert df['tax_og'].to_list() == ['4NEW1', 'arCOG0002', 'KOG0003']
    # "root" domains are replaced with the next (capitalized) level
    assert df['domain'].to_list() == ['Bacteria', 'Archaea', 'Eukaryota']


def test_parse_eggnog_annotations_polars_no_ogs(eggnog_file):
    pytest.importorskip("polars")
    df = utils.parse_eggnog_annotations_polars(eggnog_file, ogs_col=None,
                                               lazy=False)
    assert 'root_og' not in df.columns
    assert df.width == 6


SNAKE_LOG = """\
Building DAG of jobs...
Using shell: /usr/bin/bash
Job stats:
job      count
-----  -------
all          1
total        1

[Mon Jan  1 00:00:00 2024]
rule make_a:
    output: out/a.txt
    jobid: 1
    resources: tmpdir=/tmp

[Mon Jan  1 00:00:00 2024]
localrule make_b:
    output: out/b.txt
    jobid: 2
    resources: tmpdir=/tmp

[Mon Jan  1 00:00:00 2024]
checkpoint make_c:
    output: out/c.txt
    jobid: 3

[Mon Jan  1 00:00:01 2024]
Finished job 1.
1 of 4 steps (25%) done
[Mon Jan  1 00:00:02 2024]
Error in rule make_b:
    jobid: 2
    output: out/b.txt

[Mon Jan  1 00:00:00 2024]
rule all:
    input: out/a.txt, out/b.txt, out/c.txt
    jobid: 0
"""


def test_get_snake_status_table(tmp_path):
    log = tmp_path / "snakemake.log"
    log.write_text(SNAKE_LOG)

    status = utils.get_snake_status_table(str(log), as_dict=True)
    assert set(status) == {0, 1, 2, 3}
    assert status[1] == {'rule': 'make_a', 'output': 'out/a.txt', 'status': 'done'}
    assert status[2]['rule'] == 'make_b'
    assert status[2]['status'] == 'error'
    assert status[3]['rule'] == 'make_c'
    assert status[3]['status'] == 'running'
    assert status[0]['rule'] == 'all'

    table = utils.get_snake_status_table(str(log))
    assert isinstance(table, pandas.DataFrame)
    assert table.loc[1, 'status'] == 'done'
    assert set(table.columns) == {'rule', 'output', 'status'}


def test_get_snake_status_table_error_matched_by_output(tmp_path):
    """ some snakemake versions report jobid 0 for errors: match on output """
    log = tmp_path / "snakemake.log"
    log.write_text(
        "rule make_a:\n"
        "    output: out/a.txt\n"
        "    jobid: 5\n"
        "Error in rule make_a:\n"
        "    jobid: 0\n"
        "    output: out/a.txt\n"
    )
    status = utils.get_snake_status_table(str(log), as_dict=True)
    assert status[5]['status'] == 'error'


def test_get_snake_status_table_slurm(tmp_path, capsys):
    log = tmp_path / "snakemake.log"
    log.write_text(
        "rule make_a:\n"
        "    output: out/a.txt\n"
        "    jobid: 1\n"
        "Submitted job 1 with external jobid 'Submitted batch job 98765'.\n"
        "rule make_b:\n"
        "    output: out/b.txt\n"
        "    jobid: 2\n"
        "Submitted job 2 with external jobid 'Submitted batch job 98766'.\n"
        "slurmstepd: error: *** JOB 98766 ON node01 CANCELLED AT 2024-01-01 DUE TO TIME LIMIT ***\n"
        "slurmstepd: error: *** JOB 11111 ON node01 CANCELLED AT 2024-01-01 DUE TO TIME LIMIT ***\n"
    )
    status = utils.get_snake_status_table(str(log), as_dict=True, slurm=True)
    assert status[1]['slurm_id'] == 98765
    assert status[1]['status'] == 'running'
    assert status[2]['slurm_id'] == 98766
    assert status[2]['status'] == 'error'
    assert "11111" in capsys.readouterr().out


def test_get_snake_status_table_jobid_before_rule(tmp_path, capsys):
    log = tmp_path / "snakemake.log"
    log.write_text("    jobid: 1\n")
    assert utils.get_snake_status_table(str(log), as_dict=True) == {}
    assert "WARNING" in capsys.readouterr().out


###
# plotting


def test_get_N_colors(close_figures):
    colors = utils.get_N_colors(5)
    assert len(colors) == 5
    # RGBA tuples
    assert all(len(c) == 4 for c in colors)
    assert len(set(colors)) == 5

    viridis = utils.get_N_colors(3, 'viridis')
    from matplotlib import pyplot as plt
    assert viridis[0] == plt.get_cmap('viridis')(0.0)
    assert viridis[-1] == plt.get_cmap('viridis')(1.0)


def test_get_N_colors_reversed_name(close_figures):
    """ an unknown colormap name is tried spelled backwards (in reverse order) """
    assert utils.get_N_colors(4, 'sidiriv') == \
        list(reversed(utils.get_N_colors(4, 'viridis')))
    with pytest.raises(ValueError):
        utils.get_N_colors(4, 'not_a_colormap')


def test_dotplot(ax):
    hits = pandas.DataFrame({
        'qstart': [1, 501], 'qend': [500, 1000],
        'hstart': [501, 1], 'hend': [1000, 500],
        'mlen': [500, 500],
    })
    utils.dotplot(hits, ax=ax)
    lines = ax.get_lines()
    assert len(lines) == 2
    # coordinates are scaled by the longest match
    assert list(lines[0].get_xdata()) == [1 / 500, 1.0]
    assert lines[0].get_color() == 'black'

    utils.dotplot(hits, ax=ax, seq_length=False, color='red')
    assert list(ax.get_lines()[-1].get_ydata()) == [1, 500]
    assert ax.get_lines()[-1].get_color() == 'red'


def test_dotplot_makes_axes(close_figures):
    from matplotlib import pyplot as plt
    hits = pandas.DataFrame({'qstart': [1], 'qend': [10], 'hstart': [1],
                             'hend': [10], 'mlen': [10]})
    utils.dotplot(hits, subplots_kws={'figsize': [2, 2]})
    assert list(plt.gcf().get_size_inches()) == [2, 2]
    assert len(plt.gca().get_lines()) == 1


def test_save_fig_to_pdf(tmp_path, close_figures):
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    pdf_path = tmp_path / "fig.pdf"
    utils.save_fig_to_pdf(fig, str(pdf_path))
    assert pdf_path.exists()
    assert pdf_path.read_bytes().startswith(b"%PDF")
