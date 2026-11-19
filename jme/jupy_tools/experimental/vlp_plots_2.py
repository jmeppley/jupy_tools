"""
DEPRECATED: moved to jme.jupy_tools.tree_plots (TreePlotterVLP is now SyntenyTreePlotter), circular_trees and plot_metadata.

This module re-exports everything (including private names) from the new
module(s) so old notebooks keep working.
"""
import warnings as _warnings

from jme.jupy_tools import plot_metadata as _pm, circular_trees as _ct, \
    tree_plots as _tp
from jme.jupy_tools.experimental import cluster_trees, gene_synteny  # noqa
_NEW_MODULES = [_pm, _ct, _tp]

_warnings.warn(
    "jme.jupy_tools.experimental.vlp_plots_2 has moved to jme.jupy_tools.tree_plots (TreePlotterVLP is now SyntenyTreePlotter), circular_trees and plot_metadata",
    DeprecationWarning, stacklevel=2)

for _module in _NEW_MODULES:
    globals().update({_k: _v for _k, _v in vars(_module).items()
                      if not (_k.startswith('__') and _k.endswith('__'))})
del _module
