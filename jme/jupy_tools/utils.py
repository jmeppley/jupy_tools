import logging
import numpy
import pandas
import re

BLAST_COLUMNS = [
        "query",
        "hit",
        "pctid",
        "mlen",
        "mismatches",
        "gaps",
        "qstart",
        "qend",
        "hstart",
        "hend",
        "evalue",
        "score",
]
BLAST_PLUS_COLUMNS = BLAST_COLUMNS + [
    'qlen', 'hlen', 'raw_score'
]
PAF_COLUMNS = ['query','qlen','qstart','qend','strand',
               'hit', 'hlen','hstart','hend','matches',
               'mlen','mapqv']

BLAST = 'BlastTab'.upper()
BLAST_PLUS = 'BlastTab+'.upper()
PAF = 'Paf'.upper()
PAF_ALL = 'Paf+'.upper()
LASTAL = 'lastal'.upper()
LASTAL_COLUMNS = ['score',
                  'hit', 'hstart', 'hmlen', 'hstrand', 'hlen',
                  'query', 'qstart', 'qmlen', 'qstrand', 'qlen',
                  'match_string', 'eg2', 'e']


HIT_TABLE_COLUMNS = {BLAST: BLAST_COLUMNS,
                     BLAST_PLUS: BLAST_PLUS_COLUMNS,
                     PAF: PAF_COLUMNS,
                     PAF_ALL: PAF_COLUMNS + ['tp','cm','dv','rl'],
                     LASTAL: LASTAL_COLUMNS,
                    }

def computeLastHitValues(blocks):
    """
    given the query length and 'blocks' string from a last hit
    return the:
        match length

    the blocks string looks something like this:
        "73,0:1,15,0:1,13,0:1,9"
        where integer elements indicate lenghts of matches and
        colon separated elements indicate lengths of gaps
    """
    matchLen = 0
    for segment in blocks.split(','):
        try:
            matchLen += int(segment)
        except ValueError:
            (hml, qml) = segment.split(":")
            mml = max(int(hml), int(qml))
            matchLen += mml

    return matchLen


def parse_blast_m8(hit_table, format=BLAST, skiprows=0, **cutoffs):
    """ utility for quickly loading a hit table into pandas 
        cutoff keys should match column names. Use negative cutoff value to 
    """
    column_names = HIT_TABLE_COLUMNS[format]
    use_cols = list(range(len(column_names)))
    hits = \
        pandas.read_csv(
            hit_table,
            sep="\t",
            header=None,
            comment="#",
            usecols=use_cols,
            names=column_names,
            skiprows=skiprows,
        )
    
    # format specific tweaks
    if format == LASTAL:
        hits['qsmult'] = [1 if qs == '+' else -1 
                          for qs in hits.qstrand]
        hits = \
            hits.eval('hend = hstart + hmlen - 1') \
                .eval('qend = qstart + ((qmlen - 1) * qsmult)')
        hits['mlen'] = [computeLastHitValues(blocks)
                        for blocks in hits.match_string]
        hits['evalue'] = [float(re.sub('E=','',str(e)))
                          for e in hits['e']]
    if format in [PAF, PAF_ALL]:
        # calculate pctid
        hits = hits.eval('pctid = 100 * matches / mlen')
    
    query = " and ".join(f"{column} >= {value}" 
                         if value >=0 
                         else f"{column} <= {-1*value}"
                         for column, value in cutoffs.items()
                         if numpy.abs(value) > 0
                        )
    if len(query) > 0:
        return hits.query(query)
    return hits

def _print(msg, format_variables):
    if format_variables is None:
        print(msg)
    elif isinstance(format_variables, tuple):
        print(msg % format_variables)
    elif isinstance(format_variables, dict):
        print(msg.format(**format_variables))
    else:
        print(msg.format(*format_variables))


def _get_print_func(logger, level):
    if logger is None:
        return _print
    elif level <= logging.DEBUG:
        return logger.debug
    if level <= logging.INFO:
        return logger.info
    if level <= logging.WARNING:
        return logger.warning
    if level <= logging.ERROR:
        return logger.error
    return logger.critical


class LogLogger:
    """
    Class to help logging potentially long running processes
    
    It prints a message with reduced frequency in a loop so that
    the early iterations get logged, but the terminal doesn't
    get overwhelmed as the iteration progresses.
    
    If you use the old style formatting (using % as with the python logger)
    then fomratting is only performed when the message is printed.
    
        > loglogger = LogLogger(multiplier=2)
        > for x in range(int(1e6)):
        >    loglogger.log("iteration: %s", (x,))
        iteration: 0
        iteration: 1
        iteration: 3
        iteration: 7
        iteration: 15
        iteration: 31
        iteration: 63
        iteration: 127
        iteration: 255
        iteration: 511
        iteration: 1023
        iteration: 2047
        iteration: 4095
        iteration: 8191
        iteration: 16383
        iteration: 32767
        iteration: 65535
        iteration: 131071
        iteration: 262143
        iteration: 524287
    
    """
    def __init__(
        self, multiplier=2, start=0, first_trigger=1, logger=None, level=logging.INFO
    ):
        self.step = multiplier
        self.pos = start
        self.next = first_trigger
        self.__print__ = _get_print_func(logger, level)

    def check_pos(self):
        self.pos += 1
        if self.pos >= self.next:
            self.next = self.next * self.step
            return True
        return False

    def log(self, message, format_variables=None):
        """
        Reduce the frequency of log messages as the count gets higher
        """
        if self.check_pos():
            self.__print__(message, format_variables)


def get_best_fit(xd, yd, force_intercept=False, force_slope=False):
    """Return coeffs for a line of best fit"""
    # Calculate trendline
    if force_intercept:
        # intercept of 0
        x = numpy.array(xd)[:, numpy.newaxis]
        slope, _, _, _ = numpy.linalg.lstsq(x, yd, rcond=None)
        coeffs = [slope[0], 0]
        if force_slope:
            # We shouldn't get here, but let's just return the fixed values
            coeffs = (1, 0)
    elif force_slope:
        # slope of 1: intercept is average of difference
        intercept = numpy.mean(yd - xd)
        coeffs = [1, intercept]
    else:
        coeffs = numpy.polyfit(xd, yd, 1)

    return coeffs
