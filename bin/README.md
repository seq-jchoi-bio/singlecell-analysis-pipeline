# Logic Modules (bin/) — Technical README

This document summarizes the **core logic scripts** located under `bin/`:
- `IntegQC_V8.py` — Integrated scATAC‑seq Quality Control
- `FilterDoublet_V7.py` — Filtering & Doublet Removal
- `Annot_V5.py` — Cell‑type Annotation

It explains **what each module does, expected inputs, configurable globals (with defaults), outputs**, and how they interact with the top‑level wrappers (`run_qc.py`, `run_filter.py`, `run_annotation.py`). The intent is to orient users without reading the full source code.

> **Note on large artifacts:** analysis matrices (e.g., `.h5ad` / `.h5ads`) can be large and are **not included in version control**. The wrappers and logic will write them under result folders for local use.

---

## Common Design

- **Invocation:** These modules are primarily driven via the corresponding wrappers:
  - `run_qc.py` → imports and drives `bin/IntegQC_V8.py`
  - `run_filter.py` → imports and drives `bin/FilterDoublet_V7.py`
  - `run_annotation.py` → imports and drives `bin/Annot_V5.py`
- **Path injection:** Wrappers inject all paths as **plain strings** to avoid pathlib‑related edge cases in older code.
- **Genome/reference auto‑detection:** Wrappers auto‑detect `chrom.sizes`, `genes.gtf(.gz)`, and reference assets (marker CSV / atlas H5AD) under `refGenome/<species>/` (currently `rice`).
- **Result roots:** By default, wrappers write under the **current working directory**:
  - `results/QC_results/`
  - `results/Filter_results/`
  - `results/Annot_results/`
- **Entry point:** If a logic module exposes `main()`, wrappers call it; otherwise only variables are injected.

---

## `IntegQC_V8.py` — Integrated QC

### Purpose
Compute fragment‑size distributions (FSD), TSS enrichment (TSSE) metrics, and produce per‑sample QC plots and summary tables from scATAC‑seq fragment files.

### Expected inputs (set by wrapper)
- `BASE_DIR` / `base_dir`: project input root containing `<sample>/outs/fragments.tsv.gz`
- `CHROM_SIZES_PATH` / `chrom_sizes_path`: path to `*.chrom.sizes`
- `GTF_PATH` / `annot_path`: path to `genes.gtf(.gz)`
- `OUT_DIR` / `out_dir`: output directory (default `results/QC_results/`)

### Key tunables (logic defaults)
- `max_recorded_size` (default **1000**): upper bound for fragment‑size histogram
- `subsample_lines` (default **None**): subsample lines for faster preview fallback
- `grid_rows` (default **3**): rows for FSD & TSSE grid layouts
- `dpi_save` (default **300**): figure DPI

### Outputs (typical)
- `results/QC_results/fragment_size.png|svg`
- `results/QC_results/tsse_grid.png|svg`
- `results/QC_results/tsse_violin.png|svg`
- `results/QC_results/qc_summary.txt`
- `results/QC_results/qc_summary_table.csv`

---

## `FilterDoublet_V7.py` — Filtering & Doublet Removal

### Purpose
Build a minimal QC table (e.g., `is__cell_barcode`, `passed_filters=n_fragment`, `TSSE`), filter cells by TSSE and UMI bounds, and apply **Scrublet‑based doublet removal** on the retained cells. 
> **NOTE**: Scrublet 방식은 별로라고 함. 다들 검토 하길 바람.

### Expected inputs (set by wrapper)
- `base_dir`: project input root (same as above)
- `chrom_sizes_path`: path to `*.chrom.sizes`
- `annot_path`: path to `genes.gtf(.gz)`
- `out_dir`: output directory (default `results/Filter_results/`)

### Key tunables (logic defaults) 
- `TSS_CUTOFF` (default **1.0**)
- `UMI_MIN` (default **500**)
- `UMI_MAX` (default **100000**)
- `N_FEATURES` (default **250000**)
- `STRIP_BARCODE_SUFFIX` (default **False**)

### Outputs (typical)
- `results/Filter_results/merged_doublets.h5ads`  ← **large; not committed**
- `results/Filter_results/<sample>_doublets.h5ad` ← **large; not committed**
- `results/Filter_results/<sample>_filtered.csv`
- `results/Filter_results/summary.txt`

---

## `Annot_V5.py` — Cell‑type Annotation

### Purpose
Perform **marker‑based** and/or **reference‑transfer** annotation, generate UMAP visualizations, and write summary tables. Supports optional gene‑activity matrix generation.

### Expected inputs (set by wrapper)
- `BASE_DIR`: project input root
- `H5ADS_INPUT`: defaults to `BASE_DIR/Filter_results/merged_doublets.h5ads` if not overridden
- `CHROM_SIZES_PATH`: path to `*.chrom.sizes`
- `GTF_PATH`: path to `genes.gtf(.gz)`
- `REF_CSV_PATH`: marker CSV (auto‑detected under `refGenome/<species>/` if omitted)
- `REF_H5AD_PATH`: reference atlas H5AD/H5ADS (auto‑detected if omitted)
- `AGGR_DIR`: output directory (default `results/Annot_results/`)

### Key tunables (logic defaults)
- `REF_LABEL_KEY` = `"celltype_id"`
- `MARKER_MIN_LOG2FC` = **0.5**
- `MARKER_MAX_PADJ` = **0.05**
- `MIN_MARKERS_PER_SET` = **5**
- `N_FEATURES` = **250_000**
- `N_NEIGHBORS` = **10**
- `USE_DIMS` = **None** (use all spectral dims when None; may be `int` or list)
- `MAKE_GENE_ACTIVITY` = **True**
- `TRANSFER_MODE` = `"auc"` (or `"weighted"`)
- `UNKNOWN_CONF_THRESH` = **0.5**
- `SAVE_UMAP` = **True**
- `LABEL_MIN_COUNT` = **20**
- `LABEL_OUTLINE` = **True**

### Outputs (typical)
- `results/Annot_results/annotation_results.csv`
- `results/Annot_results/summary.txt`
- `results/Annot_results/Plots/umap_plots.png|svg`
- `results/Annot_results/annot_merged.h5ads`  ← **large; not committed**
- `results/Annot_results/annot_gene_activity.h5ad`  ← **large; not committed**
- `results/Annot_results/bed/*.bed` ← **large; not committed**
- `refGenome/rice/genes.gene_only.auto.gtf`

---

## Wrapper ↔ Logic Mapping (at a glance)

Wrappers only override logic defaults **when the user provides a flag**; otherwise logic defaults are preserved.

| Wrapper flag | Logic variable | Notes |
|---|---|---|
| `run_qc.py -mrs` | `max_recorded_size` | histogram range upper bound |
| `run_qc.py -sl`  | `subsample_lines`   | fast preview |
| `run_qc.py -grid`| `grid_rows`         | plot grid rows |
| `run_qc.py -dpi` | `dpi_save`          | figure DPI |
| `run_filter.py -tss` | `TSS_CUTOFF` | TSSE cutoff |
| `run_filter.py -min` | `UMI_MIN`     | |
| `run_filter.py -max` | `UMI_MAX`     | |
| `run_filter.py -feat`| `N_FEATURES`  | Scrublet/kNN features |
| `run_filter.py --strip {False/True}` | `STRIP_BARCODE_SUFFIX` | default False |
| `run_annotation.py --label` | `REF_LABEL_KEY` | reference label key |
| `run_annotation.py --log2fc` | `MARKER_MIN_LOG2FC` | |
| `run_annotation.py --pval` | `MARKER_MAX_PADJ` | |
| `run_annotation.py --markerSet` | `MIN_MARKERS_PER_SET` | |
| `run_annotation.py --feat` | `N_FEATURES` | |
| `run_annotation.py --neigh` | `N_NEIGHBORS` | |
| `run_annotation.py --dim` | `USE_DIMS` | int or comma‑list |
| `run_annotation.py --matrix {True/False}` | `MAKE_GENE_ACTIVITY` | |
| `run_annotation.py --mode {auc,weighted}` | `TRANSFER_MODE` | |
| `run_annotation.py --coreTh` | `UNKNOWN_CONF_THRESH` | |
| `run_annotation.py --umap {True/False}` | `SAVE_UMAP` | |
| `run_annotation.py --umapCount` | `LABEL_MIN_COUNT` | |
| `run_annotation.py --umapLabel {True/False}` | `LABEL_OUTLINE` | |

---

## Minimal End‑to‑End Example

```bash
# 1) QC
python run_qc.py -i /path/to/project -s rice -mrs 1500 -dpi 300

# 2) Filter + Doublet
python run_filter.py -i /path/to/project -s rice -tss 1.2 -min 800 -max 50000 --strip false

# 3) Annotation
python run_annotation.py -i /path/to/project -s rice --log2fc 0.5 --pval 0.05 --umap true
```

Outputs will be placed under `results/*` in the current working directory. Large `.h5ad`/`.h5ads` files remain locally, excluded from version control.

---

## Version & Authorship

- Logic versions: `IntegQC_V8`, `FilterDoublet_V7`, `Annot_V5`
- Wrapper version: 1.0
- Copyright (C) 2025 Sohyeong Cho, Janghyun Choi, Junbeom Lee, and Seong Kyu Han*

Please file issues or pull requests for enhancements.
