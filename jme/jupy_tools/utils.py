import logging
import io
import numpy
import pandas
import subprocess

def first(S):
    """ Simply return the first item from a collection (using next(iter(S))) """
    return next(iter(S))


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

def iterable_to_stream(iterable, str_to_bytes=True, buffer_size=io.DEFAULT_BUFFER_SIZE):
    """
    Lets you use an iterable (e.g. a generator) that yields bytestrings as a read-only
    input stream.

    The stream implements Python 3's newer I/O API (available in Python 2's io module).
    For efficiency, the stream is buffered.
    
    src: https://stackoverflow.com/a/20260030/663466
    """
    
    if str_to_bytes:
        # encode strings as bytes
        iterable = (s.encode('utf-8') for s in iterable)

    class IterStream(io.RawIOBase):
        def __init__(self):
            self.leftover = None
        def readable(self):
            return True
        def readinto(self, b):
            try:
                l = len(b)  # We're supposed to return at most this much
                chunk = self.leftover or next(iterable)
                output, self.leftover = chunk[:l], chunk[l:]
                b[:len(output)] = output
                return len(output)
            except StopIteration:
                return 0    # indicate EOF
    return io.BufferedReader(IterStream(), buffer_size=buffer_size)

def get_dataframe_from_cmd(command, shell=True, sep='\t', **kwargs):
    """ 
    Returns a pandas dataframe from a shell command that returns text in tabular format.
    
    params:
        command: the command to capture the output of
        shell: use subprocess shell mode (unsecure!) (default: True)
        sep: table delimirer (default: tab) 
        **kwargs: passed to pandas.read_csv
        
    Ideally , this would get a generator over the output lines and buffer it, but
        it's simple to just ge tthe whole output with run()"""
    
    p = subprocess.run(command, shell=shell, capture_output=True)
    return pandas.read_csv(io.BytesIO(p.stdout),
                           sep=sep, 
                           **kwargs)