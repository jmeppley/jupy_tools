"""
Recreate the layout of the Aloha_abundance_016 figure (ALOHA by depth and
month with paired size fractions, TARA by depth, metadata strips, double
color bar, log scale) from synthetic data.

    python examples/metag_heatmaps_demo.py [out.pdf|out.png]
"""
import sys

import numpy
import polars as pl

from jme.jupy_tools.experimental import metag_heatmaps
from jme.jupy_tools.utils import save_fig_to_pdf

VIRAL = metag_heatmaps.ALOHA_VIRAL        # '20-02', plotted positive
CELLULAR = metag_heatmaps.ALOHA_CELLULAR  # '00-20', plotted negative


def synthetic_data(n_genomes=30, seed=1):
    rng = numpy.random.default_rng(seed)
    genomes = [f"{'VEIME' if i % 3 else 'AFP'}_{i:03d}"
               for i in range(n_genomes)]
    depths = [25, 75, 125, 175, 250, 500, 770, 1000]
    # ALOHA: 2 filters x 8 depths x 12 months
    md = pl.DataFrame(
        [(f"S{f}_{d}_{m}", d, m, f)
         for f in (VIRAL, CELLULAR)
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
            w *= viral_bias[g] if filt == VIRAL else 1 - viral_bias[g]
            rows.append((sample, genome, scale[g] * w * rng.lognormal(0, .5)))
    hits = pl.DataFrame(rows, schema=['sample', 'genome', 'cov_per_bn_reads'],
                        orient='row')

    # TARA
    tara_samples = [f"Station{s}_{d}_0d2-3" for s in range(20, 40)
                    for d in metag_heatmaps.TARA_DEPTHS]
    tara_md = metag_heatmaps.parse_tara_samples(tara_samples)
    trows = []
    for sample, station, depth, _ in tara_md.iter_rows():
        for g, genome in enumerate(genomes):
            if rng.uniform() < .5:
                continue
            fac = {'SUR': 1, 'DCM': .6, 'MES': .2}[depth]
            if pref_depth[g] >= 500:
                fac = 1.2 - fac
            trows.append((sample, genome,
                          scale[g] * fac * rng.lognormal(0, 1)))
    tara_hits = pl.DataFrame(trows, schema=['sample', 'genome',
                                            'cov_per_bn_reads'], orient='row')
    return hits, md, tara_hits, tara_md, genomes


def main(out='metag_heatmaps_demo.png'):
    hits, md, tara_hits, tara_md, genomes = synthetic_data()

    top = metag_heatmaps.top_genomes(hits, {CELLULAR: 8, VIRAL: 15},
                                     group_col='filter', samples_md=md)
    top += [g for g in metag_heatmaps.top_genomes(tara_hits, 6)
            if g not in top]

    split = dict(split_col='filter', split_negative=CELLULAR,
                 split_order=[VIRAL, CELLULAR])
    by_depth = metag_heatmaps.aggregate(hits, md, 'depth', genomes=top,
                                        **split)
    by_month = metag_heatmaps.aggregate(hits, md, 'month', genomes=top,
                                        order=range(1, 13), **split)
    by_tara = metag_heatmaps.aggregate(
        tara_hits, tara_md, 'depth', genomes=top,
        order=list(metag_heatmaps.TARA_DEPTHS))

    (by_depth, by_month, by_tara), zero_val = metag_heatmaps.signed_log(
        [by_depth, by_month, by_tara])

    type_colors = {'AFP': 'cornflowerblue', 'VEI': 'blue'}
    type_color_fn = metag_heatmaps.color_fn_from_dict(
        type_colors, key_fn=lambda g: g[:3])
    hosts = {g: h for g, h in zip(genomes, ['Pelagibacter', None,
                                            'Cyanobacteria', 'Oceanicaulis',
                                            None] * 10) if h}
    host_colors = {'Pelagibacter': '#ffc0cb', 'Cyanobacteria': '#00ffff',
                   'Other Alphaproteobacteria': '#800080'}
    host_color_fn = metag_heatmaps.color_fn_from_dict(
        host_colors, key_fn=hosts.get,
        translate={'Oceanicaulis': 'Other Alphaproteobacteria'},
        other='#808080')

    hm = metag_heatmaps.abundance_heatmap(
        [metag_heatmaps.Panel(by_depth, 'Depth (m)\nALOHA'),
         metag_heatmaps.Panel(by_month, 'Month\nALOHA',
                              row_labels=metag_heatmaps.MONTHS),
         metag_heatmaps.Panel(
             by_tara, 'Depth\nTARA',
             row_labels=list(metag_heatmaps.TARA_DEPTHS.values()))],
        strips=[metag_heatmaps.MetadataStrip('Sequence Type', type_color_fn),
                metag_heatmaps.MetadataStrip('Mapped Host', host_color_fn)],
        label_color_fn=type_color_fn,
        log_zero_val=zero_val,
    )
    if out.endswith('.pdf'):
        save_fig_to_pdf(hm.fig, out, close=False)
    else:
        hm.fig.savefig(out, bbox_inches='tight', dpi=100)
    return hm


if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    main(*sys.argv[1:])
