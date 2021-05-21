VERSION = '0.0.1a'

# This lets me setup my env with just: from jme.jupy_tools import * 

from . import conda, hit_tables
from .filesystem import glob_wildcards
from .utils import LogLogger, first, read_tsv, get_seq_lens

from collections import Counter, defaultdict
from itertools import islice

import pandas, numpy, os, sys, re
from Bio import SeqIO
from matplotlib import pyplot as plt
import seaborn as sns
