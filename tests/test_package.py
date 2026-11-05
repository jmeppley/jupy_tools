""" package level checks and the small modules """
import re

import pandas

from jme.jupy_tools import draw_genes
from jme.jupy_tools.__version__ import VERSION


def test_version():
    # eg 0.1.5 or 0.1.5b
    assert re.match(r'^\d+\.\d+\.\d+[a-z0-9]*$', VERSION)


def test_draw_genes_is_a_stub():
    """ draw_genes is only documented so far: it accepts its args and does nothing """
    genes = pandas.DataFrame({'start': [1, 200], 'end': [150, 100],
                              'strand': ['+', '-']})
    assert draw_genes.draw_genes(genes) is None
    assert draw_genes.draw_genes(genes, location=2, thickness=.5,
                                 transform=lambda x: x,
                                 orientation='vertical') is None
