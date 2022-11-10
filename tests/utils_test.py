from jme.jupy_tools import utils

def test_log_logger():
    """ I'm not sure how to watch the logger output, but we can at least init it """
    llogger = utils.LogLogger(3, 1, 2)
    for i in range(10):
        llogger.log(f"iteration: {i}")
        
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
