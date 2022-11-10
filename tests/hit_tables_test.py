from jme.jupy_tools import hit_tables

def test_parse_m8():
    """
    def parse_blast_m8(hit_table, format=BLAST, skiprows=0, **cutoffs):
    """
    import io
    temp_file = io.StringIO()
    temp_file.write('tp1:k4.min_d0.1.nn15:d004c256-f168-4f0c-9e42-03dacdc5fd30\ttp1:k4.min_d0.5.nn15:1a8c16f8-acb3-452e-83c1-4eeab66f5a85\t96.32\t11162\t358\t28\t11150\t1\t23325\t34445\t0\t1.6e+04\t34135\t34446\t10144\ntp1:k4.min_d0.1.nn15:d004c256-f168-4f0c-9e42-03dacdc5fd30\ttp1:k5.min_d0.5.nn15:1a8c16f8-acb3-452e-83c1-4eeab66f5a85\t96.31\t11162\t358\t29\t11150\t1\t23369\t34488\t0\t1.6e+04\t34135\t34491\t10135\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp1:k6.min_d0.nn20:335cb453-9678-48d9-a9df-8da26bc90aef\t98.73\t20901\t171\t64\t34928\t14078\t1\t20857\t0\t3.15e+04\t34944\t35130\t19923\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp1:k4.min_d0.5.nn10:ec5aceaf-f0c7-41cb-bdbd-82eb386082db\t98.72\t20901\t174\t64\t34928\t14078\t1\t20857\t0\t3.15e+04\t34944\t35118\t19917\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp1:k6.min_d0.5.nn5:ec5aceaf-f0c7-41cb-bdbd-82eb386082db\t98.68\t20900\t181\t63\t34928\t14078\t1\t20855\t0\t3.15e+04\t34944\t35096\t19909\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp2:k6.min_d0.nn5:4d6b63f3-c2f3-4d55-a5b1-06194b58ed44\t98.27\t20933\t236\t72\t34928\t14078\t1\t20889\t0\t3.12e+04\t34944\t35202\t19705\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp2:k6.min_d0.nn10:4d6b63f3-c2f3-4d55-a5b1-06194b58ed44\t98.29\t20936\t228\t75\t34928\t14078\t1\t20892\t0\t3.11e+04\t34944\t35196\t19697\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp2:k5.min_d0.1.nn5:4d6b63f3-c2f3-4d55-a5b1-06194b58ed44\t98.14\t20934\t263\t69\t34932\t14078\t1\t20887\t0\t3.11e+04\t34944\t35190\t19673\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp1:k6.min_d0.nn20:335cb453-9678-48d9-a9df-8da26bc90aef\t96.31\t10323\t326\t40\t11203\t901\t24000\t34287\t0\t1.47e+04\t34944\t35130\t9281\ntp1:k4.min_d0.1.nn15:3725ea2c-3709-49f4-b879-7adf648d138a\ttp2:k6.min_d0.nn10:4d6b63f3-c2f3-4d55-a5b1-06194b58ed44\t96.56\t7961\t241\t29\t11203\t3259\t24032\t31975\t0\t1.14e+04\t34944\t35196\t7210\n')
    temp_file.seek(0)
    df = hit_tables.parse_blast_m8(temp_file, format=hit_tables.BLAST_PLUS)
    assert df.shape == (10, 15)
    assert len(set(df.hit)) == 8
    assert df.iloc[0,12] == 34135
    assert df.sort_values('pctid').pctid[0] == 96.32

def test_parse_paf():
    """
    def parse_blast_paf(hit_table, format=BLAST, skiprows=0, **cutoffs):
    """
    import io
    import numpy
    temp_file = io.StringIO()
    temp_file.write('437a607c-572a-4151-aced-2ecbaaaf2dd6\t80699\t781\t53821\t-\tcontig_18\t127941\t70093\t124878\t18878\t55497\t60\ttp:A:P\tcm:i:2404\ts1:i:18410\ts2:i:12981\tdv:f:0.0946\trl:i:747\n437a607c-572a-4151-aced-2ecbaaaf2dd6\t80699\t65401\t80609\t+\tcontig_326\t41785\t19929\t35247\t6566\t15929\t52\ttp:A:P\tcm:i:835\ts1:i:6373\ts2:i:5385\tdv:f:0.0792\trl:i:747\ne71fd889-4d97-4e3e-adad-52dbf7427f58\t75237\t55751\t75143\t+\tcontig_152\t117520\t12565\t32075\t8608\t19859\t19\ttp:A:P\tcm:i:1146\ts1:i:8492\ts2:i:8008\tdv:f:0.0760\trl:i:2049\ne71fd889-4d97-4e3e-adad-52dbf7427f58\t75237\t55679\t73934\t+\tcontig_1150\t56894\t9894\t28372\t8157\t18893\t0\ttp:A:S\tcm:i:1081\ts1:i:8008\tdv:f:0.0757\trl:i:2049\n7c41d647-ebf7-4dde-80e6-cfd18855ea93\t72268\t32537\t72194\t-\tcontig_2\t215798\t102920\t143071\t28961\t40578\t60\ttp:A:P\tcm:i:4176\ts1:i:28861\ts2:i:3414\tdv:f:0.0383\trl:i:275\n7c41d647-ebf7-4dde-80e6-cfd18855ea93\t72268\t14\t32485\t-\tcontig_577\t54763\t22091\t54755\t22826\t33066\t60\ttp:A:P\tcm:i:3302\ts1:i:22757\ts2:i:369\tdv:f:0.0409\trl:i:275\nd5b097c7-3418-4a47-b6d1-89abe4685b92\t70279\t50616\t69419\t+\tcontig_2284\t22779\t16\t19816\t7162\t19921\t0\ttp:A:P\tcm:i:932\ts1:i:6981\ts2:i:6968\tdv:f:0.0880\trl:i:201\nd5b097c7-3418-4a47-b6d1-89abe4685b92\t70279\t47390\t67102\t+\tcontig_854\t61810\t35692\t55802\t7057\t20250\t0\ttp:A:S\tcm:i:902\ts1:i:6968\tdv:f:0.0934\trl:i:201\n7e99fca7-47f8-4f8c-b060-ab3f38115fd4\t65629\t11877\t35094\t-\tcontig_1013\t89882\t36694\t60618\t2601\t24321\t60\ttp:A:P\tcm:i:289\ts1:i:2357\ts2:i:1056\tdv:f:0.1792\trl:i:2853\n7e99fca7-47f8-4f8c-b060-ab3f38115fd4\t65629\t42349\t65584\t-\tcontig_875\t161459\t100414\t124088\t2495\t24100\t60\ttp:A:P\tcm:i:281\ts1:i:2277\ts2:i:1300\tdv:f:0.1812\trl:i:2853\n5eb578b0-94df-436d-8b15-1c530394214c\t64108\t4689\t53038\t+\tcontig_1253\t108969\t59766\t108945\t23193\t49677\t60\ttp:A:P\tcm:i:3053\ts1:i:22990\ts2:i:14250\tdv:f:0.0720\trl:i:430\n5eb578b0-94df-436d-8b15-1c530394214c\t64108\t34656\t63995\t+\tcontig_341\t38096\t6187\t35393\t6081\t30031\t60\ttp:A:P\tcm:i:733\ts1:i:5840\ts2:i:1761\tdv:f:0.1332\trl:i:430\nb57d5447-df17-4619-abcd-48ea454f1eea\t63531\t3715\t52749\t+\tcontig_468\t121268\t63518\t114205\t14813\t51313\t60\ttp:A:P\tcm:i:1799\ts1:i:14354\ts2:i:6413\tdv:f:0.1063\trl:i:3871\nb57d5447-df17-4619-abcd-48ea454f1eea\t63531\t36915\t63424\t+\tcontig_557\t28416\t1619\t28327\t2151\t27315\t60\ttp:A:P\tcm:i:232\ts1:i:1908\ts2:i:1191\tdv:f:0.2010\trl:i:3871')
    temp_file.seek(0)
    df = hit_tables.parse_blast_m8(temp_file, format=hit_tables.PAF)
    assert df.shape == (14, 13)
    assert len(set(df.hit)) == 14
    assert len(set(df['query'])) == 7
    assert df.iloc[0,10] == 55497
    assert numpy.round(df.sort_values('pctid').pctid[0], 3) == 34.016
    
    
    

