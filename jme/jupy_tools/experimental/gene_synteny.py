"""
DEPRECATED: moved to jme.jupy_tools.gene_maps.

This module re-exports everything (including private names) from the new
module(s) so old notebooks keep working.
"""
import warnings as _warnings

from jme.jupy_tools import gene_maps as _new
_NEW_MODULES = [_new]

_warnings.warn(
    "jme.jupy_tools.experimental.gene_synteny has moved to jme.jupy_tools.gene_maps",
    DeprecationWarning, stacklevel=2)

for _module in _NEW_MODULES:
    globals().update({_k: _v for _k, _v in vars(_module).items()
                      if not (_k.startswith('__') and _k.endswith('__'))})
del _module
