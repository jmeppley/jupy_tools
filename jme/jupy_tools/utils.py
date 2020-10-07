import logging
import numpy
import pandas

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
HIT_TABLE_COLUMNS = {BLAST: BLAST_COLUMNS,
                     BLAST_PLUS: BLAST_PLUS_COLUMNS,
                     PAF: PAF_COLUMNS,
                     PAF_ALL: PAF_COLUMNS + ['tp','cm','dv','rl']
                    }

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
