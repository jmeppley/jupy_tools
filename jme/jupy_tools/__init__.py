# This lets me setup my env with just: from jme.jupy_tools import * 

from jme.jupy_tools import conda, hit_tables
from jme.jupy_tools.filesystem import glob_wildcards, go_to_dir
from jme.jupy_tools.utils import LogLogger, first, read_tsv, get_seq_lens, \
                                 save_fig_to_pdf

from collections import Counter, defaultdict
from itertools import islice, combinations, chain

import pandas, numpy, os, sys, re, yaml
from Bio import SeqIO
from matplotlib import pyplot as plt
import seaborn as sns