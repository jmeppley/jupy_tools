# This lets me setup my env with just: from jme.jupy_tools import * 

# the things I use most from this module
from jme.jupy_tools import conda, hit_tables
from jme.jupy_tools.filesystem import glob_wildcards, go_to_dir
from jme.jupy_tools.utils import LogLogger, first, read_tsv, \
                                 save_fig_to_pdf

# my favorite things that come with python
from collections import Counter, defaultdict
from itertools import islice, combinations, chain

import json
import os
import re
import sys
import warnings 

# these are basic enough to be required
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

# I like the following modules, but they are not always present
# let's not make them strict requirements
try:
    # I like saving things as YAML, but that's just me
    import yaml
except:
    warnings.warn(" yaml import failed. Skipping")
try:
    # super useful, but only if you are doing bioinformatics
    from Bio import SeqIO
    from jme.jupy_tools.biopython import get_seq_lens
except:
    warnings.warn(" biopython import failed. Skipping")
try:
    # a nice plitting library, but not everyone uses it
    import seaborn as sns
except:
    warnings.warn(" seaborn import failed. Skipping")

