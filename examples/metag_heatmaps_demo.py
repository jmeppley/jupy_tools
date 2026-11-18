"""
Recreate the layout of the Aloha_abundance_016 figure (ALOHA by depth and
month with paired size fractions, TARA by depth, metadata strips, double
color bar, log scale) from synthetic data.

    python examples/metag_heatmaps_demo.py [out.pdf|out.png]
"""
import sys

import numpy
import polars as pl

from jme.jupy_tools.experimental import metag_heatmaps as mh


def synthetic_data(n_genomes=30, seed=1):
    rng = numpy.random.default_rng(seed)
    genomes = [f"{'VEIME' if i % 3 else 'AFP'}_{i:03d}" for i in range(n_genomes)]
    depths = [25, 75, 125, 175, 250, 500, 770, 1000]
    # ALOHA: 2 filters x 8 depths x 12 months
    md = pl.DataFrame(
        [(f"S{f}_{d}_{m}", d, m, f)
         for f in (mh.ALOHA_VIRAL, mh.ALOHA_CELLULAR)
         for d in depths for m in range(1, 13)],
        schema=['sample', 'depth', 'month', 'filter'], orient='row')
    # each genome has a preferred depth, season and fraction
    pref_depth = rng.choice(depths, n_genomes)
    pref_month = rng.integers(1, 13, n_genomes)
    viral_bias = rng.uniform(0, 1, n_genomes)
    scale = rng.lognormal(2, 1.2, n_genomes)
    rows = []
    for sample, depth, month, filt in md.iter_rows():
        for g, genome in enumerate(genomes):
            if rng.uniform() < .35:
                continue  # absent from this sample
            w = numpy.exp(-abs(numpy.log(depth / pref_depth[g]))) * \
                (1 + .5 * numpy.cos((month - pref_month[g]) * numpy.pi / 6))
            w *= viral_bias[g] if filt == mh.ALOHA_VIRAL else 1 - viral_bias[g]
            rows.append((sample, genome, scale[g] * w * rng.lognormal(0, .5)))
    hits = pl.DataFrame(rows, schema=['sample', 'genome', 'cov_per_bn_reads'],
                        orient='row')

    # TARA
    tara_samples = [f"Station{s}_{d}_0d2-3" for s in range(20, 40)
                    for d in mh.TARA_DEPTHS]
    tara_md = mh.parse_tara_samples(tara_samples)
    trows = []
    for sample, station, depth, _ in tara_md.iter_rows():
        for g, genome in enumerate(genomes):
            if rng.uniform() < .5:
                continue
            fac = {'SUR': 1, 'DCM': .6, 'MES': .2}[depth]
            if pref_depth[g] >= 500:
                fac = 1.2 - fac
            trows.append((sample, genome, scale[g] * fac * rng.lognormal(0, 1)))
    tara_hits = pl.DataFrame(trows, schema=['sample', 'genome',
                                            'cov_per_bn_reads'], orient='row')
    return hits, md, tara_hits, tara_md, genomes


def main(out='metag_heatmaps_demo.png'):
    hits, md, tara_hits, tara_md, genomes = synthetic_data()

    top = mh.top_genomes(hits, {mh.ALOHA_CELLULAR: 8, mh.ALOHA_VIRAL: 15},
                         group_col='filter', samples_md=md)
    top += [g for g in mh.top_genomes(tara_hits, 6) if g not in top]

    split = dict(split_col='filter', split_negative=mh.ALOHA_CELLULAR,
                 split_order=[mh.ALOHA_VIRAL, mh.ALOHA_CELLULAR])
    by_depth = mh.aggregate(hits, md, 'depth', genomes=top, **split)
    by_month = mh.aggregate(hits, md, 'month', genomes=top,
                            order=range(1, 13), **split)
    by_tara = mh.aggregate(tara_hits, tara_md, 'depth', genomes=top,
                           order=list(mh.TARA_DEPTHS))

    (by_depth, by_month, by_tara), zero_val = mh.signed_log(
        [by_depth, by_month, by_tara])

    type_colors = {'AFP': 'cornflowerblue', 'VEI': 'blue'}
    hosts = {g: h for g, h in zip(genomes, ['Pelagibacter', None,
                                            'Cyanobacteria', 'Oceanicaulis',
                                            None] * 10) if h}
    host_colors = {'Pelagibacter': '#ffc0cb', 'Cyanobacteria': '#00ffff',
                   'Other Alphaproteobacteria': '#800080'}
    hm = mh.abundance_heatmap(
        [mh.Panel(by_depth, 'Depth (m)\nALOHA'),
         mh.Panel(by_month, 'Month\nALOHA', row_labels=mh.MONTHS),
         mh.Panel(by_tara, 'Depth\nTARA',
                  row_labels=list(mh.TARA_DEPTHS.values()))],
        strips=[
            mh.MetadataStrip('Sequence Type', mh.color_fn_from_dict(
                type_colors, key_fn=lambda g: g[:3])),
            mh.MetadataStrip('Mapped Host', mh.color_fn_from_dict(
                host_colors, key_fn=hosts.get,
                translate={'Oceanicaulis': 'Other Alphaproteobacteria'},
                other='#808080')),
        ],
        label_color_fn=mh.color_fn_from_dict(type_colors,
                                             key_fn=lambda g: g[:3]),
        log_zero_val=zero_val,
    )
    if out.endswith('.pdf'):
        mh.save_pdf(hm.fig, out)
    else:
        hm.fig.savefig(out, bbox_inches='tight', dpi=100)
    return hm


if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    main(*sys.argv[1:])
