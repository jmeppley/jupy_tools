# jupy_tools
A collection of assorted code that I keep coming back to. 

## Description

If I've used a block of code on more than one or two projects, I try to move it here instead of just copying it because it's likely I'll want it again in the future. This also let's me add some tests to make sure it's doing what I expect. 

This started as just a subfolder in my notebooks, but I'm spinning it out into a new project. Functions that have gotten a lot of use will be in top level modules.

The experimental folder will have newer code that's only been used a couple times and will proably change a lot as I generalize and/or improve it.

## Caveat
This repo is primarily intended for me, but I'm happy to share it. I find some of these really helpful, and maybe you will, too. However, there is no guarantee a consistent API. Stuff in experimental/ is almost guaranteed to change. Everything else *should* be more stable, but I make no promises.  

If you find anything useful, let me know and I'll either make a note to myself to try to limit API changes or I can spin it off as a stand-alone thing.

## Highlights:
Some of the things that I think others might find useful are:

### glob_wildcards
I sort of reverse-engineered the super useful glob_wildcards method from snakemake with a few differences:

 * Instead of just the wildcard values, it returns a 2-item tuple of the matching file path, and the wildcard values for that file
 * the wildcards can be returned as a dict or a NamedTuple
 * wildcards can't span multiple levels of a filesystem hierarchy
 * constraints are passed as a dict
 
### conda
Something I've struggled with in jupyter is that, while I use any conda env as a kernel for any notebook, it does't pick up the BASH execution path. So I can't. for example, just run minimap2 with a bang (!minimap2). This module does two things:

 * by just importing it, the bash execution PATH (os.environ['PATH']) gets updated to match the kernel's conda env. You will now be able to run any command line programs installed in your conda env.
 * you can 'activate' any other conda env. This is probably reckless as it alters both sys.path and environ['PATH'], but it can be useful. Only use this to run command line programs. Attempting to load python modules from other environments is almost certain to cause problems.

### nb_searcher

A tool for finding notebooks by pattern matching in file names and file
contents. It's still crude, but I find it useful.

### hit_tables
 
Tools for parsing tables in paf or blast output formats into pandas dataframes.

### imports
Importing names from the root module:

    > from jme.jupy_tools import *

gives you the following in your namespace:

```(python)
# the most useful things from this package
from jme.jupy_tools import conda, hit_tables
from jme.jupy_tools.filesystem import glob_wildcards, go_to_dir
from jme.jupy_tools.utils import LogLogger, first, read_tsv, get_seq_lens, \
                                 save_fig_to_pdf

# things I use often in jupyter
from collections import Counter, defaultdict
from itertools import islice, combinations, chain

import pandas, numpy, os, sys, re, yaml
from Bio import SeqIO
from matplotlib import pyplot as plt
import seaborn as sns
```

I dont know why I use the full names of pandas an numpy but I abbreviate
seaborn, but I do.

## Testing

I'm writing this for future me. I am testing different testing protocols with
my personal repositories, so I may not use the same approach anywhere else.

### Running nose
This is really simple

First build and activate the testing environment defined in tests/test.yml

    $ mamba env create -p ./test.env -f tests/test.yml
    $ conda activate ./test.env

Then run nosetests

    $ nosetests

### Making tests

See the `nose` documentation, but the short version is that the `nosestests`
command finds every `*_test.py` file and runs any functions starting with
`test`.
