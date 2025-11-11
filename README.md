# Single-cell ATAC-seq Integrated Pipeline
*README prepared by Janghyun Choi*

This pipeline processes `fragments.tsv.gz` files generated from Cell Ranger (scATAC-seq or scRNA-seq) pipelines, performing sequential quality control, filtering, clustering, and cell-type annotation. It is composed of four **independent wrappers (`run_*.py`)** and their corresponding **logic modules (`bin/...`)**.
Wrappers manage CLI arguments, parameter injection, and environment setup, while logic modules execute the actual analysis steps.


![Python](https://img.shields.io/badge/Python-3.10-blue?logo=python) ![micromamba](https://img.shields.io/badge/micromamba-env-green?logo=anaconda) ![conda-forge](https://img.shields.io/badge/channel-conda--forge-orange?logo=conda-forge) ![snapATAC2](https://img.shields.io/badge/snapATAC2-2.8.0-9cf?logo=rust) ![scanpy](https://img.shields.io/badge/scanpy-1.11.4-9cf?logo=python) ![status](https://img.shields.io/badge/status-stable-green)

---

## Project Layout
```
singlecell-analysis-pipeline/
├─ run_qc.py                    # Wrapper for QC
├─ run_filter.py                # Wrapper for Filter
├─ run_clustering.py            # Wrapper for Clustering
├─ run_finalReport.py           # Wrapper for Annotation
├─ bin/
│  ├─ IntegQC.py                # Logic module for QC
│  ├─ makeQCReport.py           # Logic module for QC
│  ├─ FilterDoublet.py          # Logic module for Filter
│  ├─ CLS.py                    # Logic module for Clustering
│  ├─ postCLS.py                # Logic module for Clustering
│  ├─ dotplot.py                # Logic module for Clustering
│  ├─ annotate.py               # Logic module for Annotation
│  └─ FinalReport.py            # Logic module for Annotation
├─ refGenome/                   # Custom references
│  ├─ rice/
│  │  ├─ genome.chrom.sizes
│  │  └─ genes.gtf
│  ├─ hg38/
│  │  ├─ genome.chrom.sizes
│  │  └─ genes.gtf
│  └─ make_ref.sh               # Make custom reference
├─ results/
├─ scripts/
└─ utils/                       # Other codes and utils
```

---

## Installation & Configuration

#### 1) Requirements
- OS: Linux (recommended)
- Conda/Mamba: Miniconda or Mambaforge

#### 2) Environment setup (using `environment.yml`)

Place an `environment.yml` file in the project root.  

**Create and activate:**
```bash
conda env create -f environment.yml -n {env_name}
conda activate {env_name}
```

**Verify installation:**
```bash
python -c "import snapatac2 as snap; import anndata, scanpy; import matplotlib, plotly, pandas; print('Environment OK')"
```

---

## Step-by-step Manual

> Each wrapper (`run_*.py`) supports its own CLI parameters.  
> Below is a concise overview of key arguments for each step.  
> For detailed usage and advanced options, **run each script with the `-h` or `--help` flag**.

### Working Directory Recommendation

- While this pipeline allows **independent specification of input and output paths** for each step, it is **strongly recommended to execute all wrappers from the project root directory**.
- By default, each wrapper automatically creates the following subdirectories:

```     
    Root
    ├─ `QC_results/`
    ├─ `Filter_results/`
    ├─ `Annot_results/`
    └─ `FinalReport.html`
```

- These directories are referenced internally by downstream scripts (e.g., plotting, report generation, summary merging).
- Running all steps from a single root ensures that file paths remain consistent and that the final report and auxiliary analyses function correctly.

---

### Step 1 — Quality Control (QC)

**Purpose:**  
Generate fragment-size distributions (FSD), TSSE metrics, and per-sample QC grids.  
Output summary tables and an optional HTML report.

**Input:**  
- `<sample>/outs/fragments.tsv.gz` automatically detected.  
- Reference files: `chrom.sizes`, `genes.gtf` 
> **If a user’s reference genome is not available,**
this pipeline allows converting the gtf and fasta files used in Cell Ranger into a compatible format for direct use.
**To build a custom reference,** refer to `refGenome/README,md` or use the `-h` flag of `refGenome/make_ref.sh`.

**Core parameters:**
| Parameter | Description |
|------------|--------------|
| `-i`, `--input` | Base directory containing sample folders with `outs/fragments.tsv.gz`. |
| `-o`, `--outdir` | Output directory (default: `./QC_results`). |
| `-s`, `--species` | Species ID for locating `refGenome/<species>/genes.gtf` and `genome.chrom.sizes`. |
| `-j`, `--n-jobs` | Number of parallel threads. |

> See `python run_qc.py -h` for the complete argument list.

**Output:**  
- `QC_results/` with PNG/SVG plots, text summaries, and optional HTML report.

**Example:**  
```bash
python run_qc.py -i ./data -s hg38 -j 10
```

---

## Step 2 — Filtering and Doublet Removal

**Purpose:**  
Perform TSSE/UMI-based filtering and Scrublet doublet detection.

**Input:**  
- `<sample>/outs/fragments.tsv.gz` automatically detected.  
- Reference files: `chrom.sizes`, `genes.gtf`

**Core parameters:**
| Parameter | Description |
|------------|--------------|
| `-i`, `--input` | Input directory or h5ad file from Step 1. |
| `-s`, `--species` | Reference species (for GTF/chrom.sizes lookup) from Step 1. |
| `--tss` | Minimum TSS enrichment threshold. |
| `--min`, `--max` | Min/Max fragment counts per cell. |
| `-j`, `--n-jobs` | Threads for filtering and Scrublet. |

> See `python run_filter.py -h` for details.

**Output:**  
- `Filter_results/` folder with filtered data (`*.h5ad/h5ads`) and `summary.txt`.

**Example:**  
```bash
python run_filter.py -i ./data -s hg38 -tss 10.0 -min 5000 -max 100000
```

---

## Step 3 — Clustering

**Purpose:**  
Conduct feature selection, spectral embedding, Harmony integration, Leiden clustering, and UMAP visualization.

**IMPORTANT!!**
- This pipeline does not rely on the automatic `snap.tl.make_gene_matrix()` function of SnapATAC.
- Instead, it directly processes the input reference genome (GTF and chrom.sizes) to build the gene-activity matrix,
enabling compatibility with custom or non-model species.

**Input:**  
- `Filter_results/merged_doublets.h5ads` automatically detected.  
- Reference files: `chrom.sizes`, `genes.gtf`
- Marker genes: A `csv` file containing marker genes. **The file must include the headers `type` and either `gene_id` or `gene_name`.**

**Core parameters:**
| Parameter | Description |
|------------|--------------|
| `-i`, `--input` | Input h5ads file (auto-detected obtained from Step 2). |
| `-s`, `--species` | Reference genome from Step 1. |
| `-r`, `--markers` | Marker gene list as a csv format. |

> See `python run_clustering.py -h` for advanced options.

**Output:**  
- All results are saved under the automatically generated directory `Annot_results/` in the current working folder.
- `Annot_results/plots/`: UMAP plots and coordinate data.
- `Annot_results/clustered.h5ad`: Cluster-level feature matrix.
- `Annot_results/annot_gene_activity.h5ad`: Matrix combining clustering and Leiden information.
- `Annot_results/dotplots/`: Dotplot figures and numerical data tables.

**Example:**  
```bash
python run_clustering.py -i ./ -s hg38 -r humanMarker.csv
```

---

## Step 4 — Annotation and Final Report

**Purpose:**  
Perform marker-based annotation, update the annotated h5ad file with cell-type information, and export BED files per cell type, followed by generating a unified final report that integrates the dotplot, label map, and UMAP.

**Input:**  
- `QC_results/summary.txt`  
- `Filter_results/summary.txt`  
- `Annot_results/plots/umap.csv`  
- `Annot_results/dotplots/*.csv`  
- **Key input**: **`Annot_results/dotplots/label_map.csv`**
> - This file plays a critical role in the annotation process.
> - During **Step 3**, Leiden clusters are automatically mapped to cell types based on marker-pattern signatures a method conceptually similar to `scanpy.tl.score_genes`.
> - The resulting label_map.csv file contains two columns: leiden and type.
> - **If the automatic assignment does not match the user’s interpretation,**
the type column can be manually edited to redefine the cell-type labels for each Leiden cluster.
> - **Step 4** directly uses this updated label_map.csv to perform annotation,
making this file the **central bridge between automated scoring and user-defined curation** in the pipeline.

**Core parameters:**
| Parameter | Description |
|------------|--------------|
| `-i`, `--input` | Base root project directory (with `Annot_results/` and `Filter_results/`). |
| `-r`, `--marker` | Optional marker gene CSV for generating report. |
| `-n`, `--name` | Custom report name (optional, default: `FinalReport.html`). |

> See `python run_AnnotReport.py -h` for CLI reference.

**Output:**  
- `Annot_results/update_h5ad/`: Update the filtered h5ad file with cell-type information/
- `Annot_results/bed/`: Export bed files per cell-type.
- `./FinalReport.html`: Integrated HTML report (full embeded).

**Example:**  
```bash
python run_AnnotReport.py -i ./ -r humanMarker.csv -n HumanData.html
```

---

## Version & Authorship

| Step | Wrapper |  Logic Module(s) | Version |
|------|----------|----------------|----------|
| **1 – QC process** | `run_qc.py` | `bin/IntegQC.py`, `bin/makeQCReport.py` | v2.2 |
| **2 – Filter/Doublet** | `run_filter.py` | `bin/FilterDoublet.py` | v3.0 |
| **3 – Clustering** | `run_clustering.py` | `bin/CLS.py`, `bin/postCLS.py`, `bin/dotplot.py` | v1.3 |
| **4 – Annotation/Report** | `run_AnnotReport.py` | `bin/annotate.py`, `bin/FinalReport.py` | v4.0 |

**Main Developer and Code Architect:**  
- Janghyun Choi, Ph.D. — implemented the main pipeline modules and designed major functional components.

**Project Supervisor:**
- Seong Kyu Han, Ph. D. — provided conceptual guidance and overall project supervision.

**Contributors:**  
- Sohyeong Cho — contributed to code implementation, module development, and pipeline testing.  
- Junbeom Lee – supported data processing, logic refinement, and validation.
