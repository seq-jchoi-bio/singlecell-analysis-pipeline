# Single-cell ATAC-seq Integrated Pipeline
**[[github guide](README_github.md)]**, **[[Logic code manual](bin/README.md)]**
This repository provides three **command-line wrappers** that serve as main entry points for the single-cell ATAC-seq pipeline.
Each wrapper reads user arguments, auto-detects input format and reference genomes (Currently, only rice genomes are supported.), and passes configuration to the corresponding logic scripts inside `bin/`.

![Python](https://img.shields.io/badge/Python-3.10-blue?logo=python) ![micromamba](https://img.shields.io/badge/micromamba-env-green?logo=anaconda) ![conda-forge](https://img.shields.io/badge/channel-conda--forge-orange?logo=conda-forge) ![snapATAC2](https://img.shields.io/badge/snapATAC2-2.8.0-9cf?logo=rust) ![status](https://img.shields.io/badge/status-developing-yellow)

## Project Layout
```
singlecell-analysis-pipeline/
├─ run_qc.py          # Quality Control (Integrated scATAC-seq QC)
├─ run_filter.py      # Filter + Doublet Removal
├─ run_annotation.py  # Cell-type Annotation
└─ bin/               # Logic codes
   ├─ IntegQC_V8.py
   ├─ FilterDoublet_V7.py
   └─ Annot_V5.py
└─ refGenome/        # Custom reference
   └─ rice/
       ├─ chrom.sizes
       ├─ genes.gtf
       └─ sigle-cell_marker.csv
└─ results/
   ├─ QC_results/             # run_qc.py
   ├─ Filter_results/         # run_filter.py
   └─ Annotation_results/     # run_annotation.py
└─ utils/                     # Other codes and utils
```

## Output Overview
Each process stores its results in the `results/` folder.  

| Process | Output folder | Main files (examples) | Description |
|---------|--------------|------------------------|------------|
| **QC** | `results/QC_results/` | `fragment_size.png/svg`, `TSSE_violin.png/svg`, `TSSE_grid.png/svg`, `qc_summary(_table).txt/csv` | Fragment-size distributions, TSSE metrics, combined QC plots, and text/csv summary |
| **Filter + Doublet** | `results/Filter_results/` | `summary.txt`, `<sample>_filtered.csv`, `<sample>_doublets.h5ad`, `merged_doublets.h5ads` | QC tables, filtered cell barcodes, merged h5ad files |
| **Annotation** | `results/Annotation_results/` | `Plots/*.png/svg`, `bed/*.bed`, `annot_gene_activity.h5ad`, `annot_merged_cells.h5ad`, `annotation_results.csv`, `summary.txt` | cluster annotations, UMAP visualization |

> **Note**: `.h5ad` and other large files generated during analysis are automatically placed under these result folders for internal use but are **excluded from the repository** to reduce storage size.

## Example File
The example data used in the code are located in the server at **`/data/pipelines/snatac/Example`**.
This dataset contains a chromosome 1 fragment file from scATAC-seq data, organized as follows:

```
/data/pipelines/snatac/Example/
└── Sample1/
    └── outs/
        └── fragments.tsv.gz
└── Sample2/
    └── outs/
        └── fragments.tsv.gz
```

## Step 1. Quality Control
Run QC as follows:
```bash
python run_qc.py -i <INPUT_DIR> -s rice
```
Main options:
- `-mrs {NUMBER}` : Upper bound for fragment-size histogram (default: 1000)
- `-sl {NUMBER}` : Subsample number of lines for fast preview
- `-grid {NUMBER}` : Grid rows for plots (default: 3)
- `-dpi {NUMBER}` : Figure DPI (default: 300)

## Step 2. Filter + Doublet Removal
Run filtering and doublet detection:
```bash
python run_filter.py -i <INPUT_DIR> -s rice -tss <float> -min <int> -max <int>
```
> **Note**: Marker genes currently support only CSV and h5ad file formats and are located in the `refGenome/`.

Key options:
- `-tss {FLOAT}` : TSSE cutoff (default: 1.0)
- `-min {INT}` : Minimum UMI (default: 500)
- `-max {INT}` : Maximum UMI (default: 100000)
- `-feat {INT}` : Number of features for scrublet/kNN graph (default: 250000)
- `--strip {False/True}` : Strip trailing digits from barcodes (default: False)

## Step 3. Cell-type Annotation
Run cell-type annotation with reference-based transfer:
```bash
python run_annotation.py -i <INPUT_DIR> -s rice
```
Key options:
- `--label {KEY}` : Map to REF_LABEL_KEY (default: celltype_id)
- `--log2fc {FLOAT}` : Min log2 fold-change for marker genes (default: 0.5)
- `--pval {FLOAT}` : Max adjusted p-value (default: 0.05)
- `--markerSet {INT}` : Min markers per set (default: 5)
- `--feat {INT}` : Number of features (default: 250000)
- `--neigh {INT}` : Number of neighbors (default: 10)
- `--dim {INT/None}` : Spectral dimensions to use (default: None)
- `--matrix {True/False}` : Generate gene activity matrix (default: True)
- `--mode {auc/weighted}` : Transfer mode (default: auc)
- `--coreTh {FLOAT}` : Confidence threshold (default: 0.5)
- `--umap {True/False}` : Save UMAP plot (default: True)
- `--umapCount {INT}` : Min count per UMAP label (default: 20)
- `--umapLabel {True/False}` : Outline labels on UMAP (default: True)

## General Notes
- All wrappers require `-i` (input directory) and `-s` (species).
- Species auto-detection supports `rice` and can be extended to other genomes by adding folders under `refGenome/`.
- If an argument is omitted, each wrapper falls back to the default value **defined in the logic scripts** inside `bin/`.

## Version & Authorship

- Logic and wrapper versions: `IntegQC_V8 (version 2.2)`, `FilterDoublet_V7 (version 2.2)`, `Annot_V5 (version 2.8)`
- Wrapper version: 1.0
- Copyright (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*

Please file issues or pull requests for enhancements.