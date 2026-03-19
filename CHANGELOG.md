# Changelog
## 1.3.0
* Add `bean fuse` subcommand — computes FUSE (Functional Score Using Structural Ensemble) scores from a `bean run` element result CSV
  * James-Stein positional shrinkage estimator + FUNSUM substitution matrix scoring
  * Optional secondary-structure-specific FUNSUM scoring via DSSP `.dss` file (`--dss`)
  * `--mu-sd-max` flag to filter low-confidence variants before scoring
  * `--exclude-lof` flag to omit stop-gain variants and use the noLOF FUNSUM
  * Automatically uses `mu_z_adj` when available, falls back to `mu_z` with a warning
  * Handles single-letter and three-letter amino acid codes; parses `[gene:]pos:ref>alt` and `CodingNoncodingAllele` edit string formats
  * Bundles pre-converted FUNSUM matrices as pickles (`bean/fuse/data/`) — no R dependency at runtime (requires `pyreadr` and `rdata` for one-time conversion only)
* Add `crispr_bean_workflow.ipynb` — end-to-end tutorial notebook covering `count-samples` → `qc` → `filter` → `run` → `fuse` with plotnine visualisations
## 1.2.9
* Add sgRNA information as output
* Add pseudocounts for LFC calculation (`perturb-tools>=0.3.5`)
## 1.2.8
* Change .pyx files to be compatible with more recent numpy versions
## 1.2.7
* **CRITICAL** Fix sample ordering & masking issue for survival screens
## 1.2.6
* Fix overflow in `bean run survival` and autograde error related to inplace assignment for `bean run survival tiling`.
## 1.2.5
* Allow `bean run .. tiling` for untranslated `--allele-df-key`.