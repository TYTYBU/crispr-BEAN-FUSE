# <img src="docs/assets/bean_title2.svg" alt="crispr-bean" height="50"/>

[![PyPI pyversions](https://img.shields.io/pypi/pyversions/crispr-bean)](https://pypi.org/project/crispr-bean/)
[![PyPI version](https://img.shields.io/pypi/v/crispr-bean)](https://pypi.org/project/crispr-bean/)
[![Test](https://github.com/pinellolab/crispr-bean/actions/workflows/CI.yml/badge.svg)](https://github.com/pinellolab/crispr-bean/actions/workflows/CI.yml)
[![Documentation](https://github.com/pinellolab/crispr-bean/actions/workflows/documentation.yml/badge.svg)](https://github.com/pinellolab/crispr-bean/actions/workflows/documentation.yml)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

Variant effect estimation from base editing screens, with structure-informed score denoising via FUSE.

<img src="docs/assets/summary.png" alt="Reporter construct" width="700"/>

## Installation
First install [PyTorch](https://pytorch.org/get-started/), then:
```bash
pip install crispr-bean
```

For the latest version (and test data), install from source:
```bash
git clone https://github.com/pinellolab/crispr-bean.git
cd crispr-bean && pip install -e .
```

## bean run
Quantify variant effect sizes from screen data.

```bash
# Variant library (e.g. GWAS SNPs)
bean run variant   my_screen.h5ad --output-prefix results/

# Coding tiling library
bean run tiling    my_screen.h5ad --output-prefix results/

# Survival / proliferation screen
bean run survival  my_screen.h5ad --output-prefix results/
```

See the [model](https://pinellolab.github.io/crispr-bean/model.html) and [output format](https://pinellolab.github.io/crispr-bean/run.html#output) documentation for details.

## bean fuse
Denoise per-variant scores from `bean run` using the FUNSUM substitution matrix and optional protein secondary structure (FUSE — Functional Score Using Structural Ensemble).

```bash
# Basic (auto-detects mu_z_adj or mu_z as raw score)
bean fuse bean_element_result.MixtureNormal.csv -o results/

# With secondary-structure-specific FUNSUM and confidence filter
bean fuse bean_element_result.MixtureNormal.csv \
  --dss structure.dss --dssp-offset 21 \
  --mu-sd-max 1.0 -o results/

# If bean run was run without bean filter (target column not yet AA-translated)
bean fuse bean_element_result.MixtureNormal.csv \
  --target-decoder LDLR_target_decoder.csv \
  --dss structure.dss --dssp-offset 21 \
  -o results/
```

Output columns: `gene`, `aapos`, `aaref`, `aaalt`, `functional_class`, `raw_score`, `norm_raw_score`, `pos_score`, `sub_score`, `sub_score_ss`, `FUSE_score`, `FUSE_SS_score`.

Example files for the LDLR peDMS dataset are in [`examples/LDLR/`](examples/LDLR/).

## Full pipeline

<img src="docs/assets/dag_bean_v2.svg" alt="dag_bean_v2.svg" height="500"/>

| Step | Command | Description |
|------|---------|-------------|
| 1 | [`count` / `count-samples`](https://pinellolab.github.io/crispr-bean/count.html) | Map guides (and reporter) from FASTQ |
| 2 | [`profile`](https://pinellolab.github.io/crispr-bean/profile.html) | Profile editor preferences |
| 3 | [`qc`](https://pinellolab.github.io/crispr-bean/qc.html) | Quality control and sample/guide filtering |
| 4 | [`filter`](https://pinellolab.github.io/crispr-bean/filter.html) | Filter alleles; translates coding variants for tiling screens |
| 5 | [`run`](https://pinellolab.github.io/crispr-bean/run.html) | Estimate variant effect sizes |
| 6 | `fuse` | Denoise scores with FUNSUM + optional secondary structure |

For a complete walkthrough see [`crispr_bean_workflow.ipynb`](crispr_bean_workflow.ipynb).
For tutorials by screen type see the [documentation](https://pinellolab.github.io/crispr-bean/).

### Variant vs tiling library design

**Variant** — several gRNAs per target variant; bystander effects ignored. Best for noncoding GWAS screens.
**Tiling** — dense gRNA tiling of a region; all alleles considered and translated to coding variants. Best for coding sequence DMS-style screens.

## Python API
```python
import bean as be
cdata = be.read_h5ad("bean_counts_sample.h5ad")
```
See the [ReporterScreen API tutorial](docs/ReporterScreen_api.ipynb) for details.

## Citation
Ryu, J., Barkal, S., Yu, T. et al. Joint genotypic and phenotypic outcome modeling improves base editing variant effect quantification. *Nat Genet* (2024). https://doi.org/10.1038/s41588-024-01726-6

## Contributing
See [CHANGELOG](CHANGELOG.md) for recent updates. Questions and PRs welcome — please open an issue.
