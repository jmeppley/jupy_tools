from jme.jupy_tools import utils
import logging, os
test_log_file = '.logger.test.out'
try:
    os.remove(test_log_file)
except FileNotFoundError:
    pass

def test_log_logger():
    """ I'm not sure how to watch the logger output, but we can at least init it """
    logger = logging.getLogger('test')
    file_handler = logging.FileHandler(test_log_file)
    logger.addHandler(file_handler)
    logger.propagate = False
    logger.warning("Init logger")
    llogger = utils.LogLogger(3, 1, 2, logger=logger, level=logging.WARNING)
    for i in range(10):
        llogger.log(f"Iteration: {i}")

    with open(test_log_file) as lines:
        lines = list(lines)
    assert len(lines) == 3
    assert lines[0] == 'Init logger\n'
    assert lines[2] == 'Iteration: 4\n'
    assert lines[1] == 'Iteration: 0\n'
    os.remove(test_log_file)
        
def test_best_fit():
    """
    def get_best_fit(xd, yd, force_intercept=False, force_slope=False):
    """
    import numpy
    
    x = numpy.array([0,1,2,3,4])
    for i in range(1,10,2):
        y = i * x
        slope, intercept = utils.get_best_fit(x, y, True)
        assert abs(slope - i) < 1e-10
        assert intercept == 0
        
        slope, intercept = utils.get_best_fit(x, y)
        assert abs(slope - i) < 1e-10
        assert abs(intercept) < 1e-10 
