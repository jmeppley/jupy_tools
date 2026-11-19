# Migrating smart_labels, distances, and mapping out of experimental

Status: **parked**. This work is low priority and is kept on its own branch
(`migrate-smart-labels-distances-mapping`) to merge later.

## What changed

Three modules moved from `jme/jupy_tools/experimental/` to
`jme/jupy_tools/`. They were picked because each stands alone and isn't tied
to a single project's data. The larger, interdependent modules
(`cluster_trees`, `vlp_plots*`, `gene_synteny`, etc.) were left in
`experimental`.

The branch has two commits:

1. **Move only:** a pure `git mv` of the three files, so the history and
   the refactor diff are easy to follow.
2. **Clean up and test:** the refactors, bug fixes, tests, and README notes
   described below.

### smart_labels.py

Labels points on a matplotlib plot, moving labels that would overlap.

- **Bug fix:** `move_bounds` crashed when `axis` was `None` (the default),
  because it called `axis.lower()`.
- `axis` is validated ('x', 'y', None, or blank).
- `arrowprops` is no longer a shared mutable default
  (`DEFAULT_ARROWPROPS` holds the defaults).
- Debug output uses `logging` instead of a `debug=True` flag.
- Raises `TooManyAdjustmentsError` (with `max_moves` configurable) instead of
  printing and raising a bare `Exception`.
- Moved labels are shifted by the total move from their original text
  position, so text alignment settings are respected.
- `annotate_points` returns the final `Annotation` objects.
- `experimental/sunburst.py` already imported `jme.jupy_tools.smart_labels`,
  which didn't exist before this move.

### distances.py

Runs a pairwise function over every pair of items with multiprocessing.

- **Thread cap:** the default worker count is all available CPUs, but at most
  `MAX_DEFAULT_THREADS` (8). Before, it started `chunk_root²` processes
  (25 by default). `threads=1` runs in-process without multiprocessing.
- `chunkify_data` no longer yields empty chunks. Tests check that every pair
  is covered exactly once for a range of sizes.
- `gene_annots_to_data_array` was rewritten. The old version assigned a
  column to a slice of the input (chained assignment) and used a
  `(dict, sum)` aggregation. It now returns a list of
  `(genome, {annot: length}, total_length)` tuples.
- `calc_shared_gene_length_ratio` takes `ignore_annots` (default
  `("Unknown",)`) and returns 0 for genomes with no genes. Its docstring notes
  that it is a similarity, not a distance.
- `data_kws` no longer uses a mutable default.

### mapping.py

Turns read mappings into per-base coverage arrays.

- **Dependency removed:** it no longer uses `edl.blastm8` (py-metagenomics).
  SAM records and CIGAR strings are parsed directly
  (`parse_sam_alignments`, `get_cigar_ref_length`).
- BAM/CRAM files are read with `samtools view -h` run as an argument list
  (no shell). A non-zero exit status raises `CalledProcessError`.
- Flag filtering uses `SAMFlag`. `DEFAULT_EXCLUDE_FLAGS` is the same as
  `samtools view -F 2308`.
- New `get_deltas_from_hit_table` works on `hit_tables` DataFrames
  (including 0-based PAF/BED coordinates).
- **Bug fixes in `deltas_to_cov`:**
  - `numpy.NAN` doesn't exist in numpy 2.
  - An `x_max` shorter than the deltas crashed.
  - `portion=True` rounding could leave tiny non-zero coverage.

  It is now vectorized (`numpy.add.at` + `cumsum`).

#### API changes in mapping

- `parse_deltas(hit_table, ...)` became `get_deltas(alignments, portion=False)`,
  taking `(read, reference, start, end)` tuples.
- `parse_deltas_from_sam` takes `exclude_flags` and `use_samtools` instead of
  a samtools command template and blastm8 filter kwargs.
- Without `x_max`, `deltas_to_cov` stops at the last covered base. It used to
  add one trailing zero position.
- `get_deltas` returns plain dicts, not defaultdicts.

## Tests

New files: `tests/test_smart_labels.py`, `tests/test_distances.py`, and
`tests/test_mapping.py`. They cover 97–100% of the three modules. The samtools
test is skipped if `samtools` isn't installed.

When this was written, the full suite had 2 failures, both already failing
before this work: `test_parse_blast` and `test_parse_blast_cutoffs` in
`tests/test_hit_tables.py`. With pandas 3, an e-value of `1e-50` is parsed as
`9.999999999999999e-51` and the tests compare with `==`.

## Before merging

- Rebase onto the current default branch and re-run `pytest`.
- Search notebooks for imports of `jme.jupy_tools.experimental.mapping`,
  `.distances`, or `.smart_labels` and update them. Also update callers of
  `mapping.parse_deltas`, whose signature changed.
- Consider whether `samtools` should be added to `tests/test.yml`.
