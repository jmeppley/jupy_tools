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

def test_parse_m8():
    """
    def parse_blast_m8(hit_table, format=BLAST, skiprows=0, **cutoffs):
    """
    import io
    temp_file = io.StringIO()
    temp_file.write('tp1:k4.min_d0.1.nn15:d004c256-f168-4f0c-9e42-03dacdc5fd30\ttp1:k4.min_d0.5.nn15:1a8c16f8-acb3-452e-83c1-4eeab66f5a85\t96.32\t11162\t358\t28\t11150\t1\t23325\t34445\t0\t1.6e+04\t34135\t34446\t10144\ntp1:k4.min_d0.1.nn15:d004c256-f168-4f0c-9e42-03dacdc5fd30\ttp1:k5.min_d0.5.nn15:1a8c16f8-acb3-452e-83c1-4eeab66f5a85\t96.31\t11162\t358\t29\t11150\t1\t23369\t34488\t0\t1.6e+04\t34135\t34491\t10135\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp1:k6.min_d0.nn20:335cb453-9678-48d9-a9df-8da26bc90aef\t98.73\t20901\t171\t64\t34928\t14078\t1\t20857\t0\t3.15e+04\t34944\t35130\t19923\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp1:k4.min_d0.5.nn10:ec5aceaf-f0c7-41cb-bdbd-82eb386082db\t98.72\t20901\t174\t64\t34928\t14078\t1\t20857\t0\t3.15e+04\t34944\t35118\t19917\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp1:k6.min_d0.5.nn5:ec5aceaf-f0c7-41cb-bdbd-82eb386082db\t98.68\t20900\t181\t63\t34928\t14078\t1\t20855\t0\t3.15e+04\t34944\t35096\t19909\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp2:k6.min_d0.nn5:4d6b63f3-c2f3-4d55-a5b1-06194b58ed44\t98.27\t20933\t236\t72\t34928\t14078\t1\t20889\t0\t3.12e+04\t34944\t35202\t19705\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp2:k6.min_d0.nn10:4d6b63f3-c2f3-4d55-a5b1-06194b58ed44\t98.29\t20936\t228\t75\t34928\t14078\t1\t20892\t0\t3.11e+04\t34944\t35196\t19697\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp2:k5.min_d0.1.nn5:4d6b63f3-c2f3-4d55-a5b1-06194b58ed44\t98.14\t20934\t263\t69\t34932\t14078\t1\t20887\t0\t3.11e+04\t34944\t35190\t19673\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp1:k6.min_d0.nn20:335cb453-9678-48d9-a9df-8da26bc90aef\t96.31\t10323\t326\t40\t11203\t901\t24000\t34287\t0\t1.47e+04\t34944\t35130\t9281\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp2:k6.min_d0.nn10:4d6b63f3-c2f3-4d55-a5b1-06194b58ed44\t96.56\t7961\t241\t29\t11203\t3259\t24032\t31975\t0\t1.14e+04\t34944\t35196\t7210\n')
    temp_file.seek(0)
    df = utils.parse_blast_m8(temp_file, format=utils.BLAST_PLUS)
    assert df.shape == (10, 15)
    assert len(set(df.hit)) == 8
    assert df.iloc[0,12] == 34135
    assert df.sort_values('pctid').pctid[0] == 96.32
    
    
    
