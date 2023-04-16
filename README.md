# ER+ Breast Cancer Immune Landscape - scRNA-seq Analysis

Analysis code, figure-level results, and manuscript materials supporting the
single-cell RNA-seq sections of:

> Onkar S, Cui J, Zou J, Cardello C, Cillo AR, Uddin MR, Sagan A, Joy M,
> Osmanbeyoglu HU, Pogue-Geile K, McAuliffe PF, Lucas PC, Tseng GC, Lee AV,
> Bruno TC, Oesterreich S, Vignali DAA.
> **Immune landscape in estrogen receptor positive breast cancer reveals a
> divergent macrophage-driven microenvironment.**
> *Nature Cancer*. 2023;4(4):516-534.
> DOI: [10.1038/s43018-023-00527-w](https://doi.org/10.1038/s43018-023-00527-w) *
> PMC: [PMC11194444](https://pmc.ncbi.nlm.nih.gov/articles/PMC11194444/)

```bibtex
@article{onkar2023immune,
  title   = {Immune landscape in estrogen receptor positive breast cancer reveals a divergent macrophage-driven microenvironment},
  author  = {Onkar, Sayali and Cui, Jian and Zou, Jian and Cardello, Carly and Cillo, Anthony R and Uddin, Mostofa Rafid and Sagan, April and Joy, Marion and Osmanbeyoglu, Hatice U and Pogue-Geile, Katherine and McAuliffe, Priscilla F and Lucas, Peter C and Tseng, George C and Lee, Adrian V and Bruno, Tullia C and Oesterreich, Steffi and Vignali, Dario AA},
  journal = {Nature Cancer},
  year    = {2023},
  volume  = {4},
  number  = {4},
  pages   = {516--534},
  doi     = {10.1038/s43018-023-00527-w}
}
```

The paper's central finding: in ER+ breast cancer, macrophages - not T cells -
are the tumor-infiltrating immune subset with the greatest transcriptional
divergence between invasive ductal carcinoma (IDC) and invasive lobular
carcinoma (ILC). This repository covers the computational analysis behind
**Figures 4-6** (clustering/annotation of the CD45+ immune compartment,
monocyte/macrophage sub-clustering, RNA velocity/pseudotime, differential
expression, SPIA pathway enrichment, and macrophage-T cell receptor-ligand
communication). Flow cytometry, mIHC, and clinical/survival analyses (Figures
1-3, 7 and most Extended Data figures) were done outside this codebase and
are not included here.

Raw sequencing data (10x gene expression matrices, CITE-seq/hashing counts,
velocyto loom files) and large intermediate Seurat/CellChat `.rds`/`.RData`
objects (0.1-4 GB each) are **not** included in this repository - GitHub is
not an appropriate place for them. They remain on the lab's OneDrive project
folder, organized into `data/` (raw, single canonical copy) and `temp/`
(intermediate analysis objects needed to re-run the pipeline).

## Repository layout

```
scripts/
  pipeline/         Core R/Rmd pipeline: initialization -> QC -> clustering ->
                     cluster labeling -> Bhattacharyya distance -> monocyte/
                     macrophage sub-clustering -> differential expression ->
                     SPIA pathway enrichment -> celltalker receptor-ligand
                     analysis -> CellChat exploration
  velocity/          scVelo (Python) RNA velocity + latent-time analysis
  preprocessing/      velocyto raw-data preprocessing scripts (.sh/.sbatch)

results/
  figure4_immune_landscape_monocyte_macrophage/
  figure5_monomacro_trajectory_DE/
  figure6_pathway_cell_communication/
    each contains:
      source_data/    official per-panel source data (as deposited with the paper)
      figures/         rendered plots (png/pdf) produced by the scripts above
      *_tables/        supporting DE/GO/KEGG/receptor-ligand tables

manuscript/
  method.docx                                              methods draft (RNA velocity + CellTalker sections)
  TonysMarker List_differentially_expressed_genes.xlsx      DE marker gene reference list
  Sylie101921-Curated gene list for scRNAseq[4895].docx     curated gene list for scRNA-seq annotation
```

## Notes on scope

- `scripts/pipeline/` holds the versions of the analysis that fed the
  published figure panels (the project's `2.Analysis/` pipeline, iterated on
  from 2021-2023). An earlier 2021 draft pass of the same pipeline (root-level
  `1.test.R`, `2.clustering.R`, `4.x`, `5.x` scripts) was superseded by this
  pipeline and is intentionally excluded, along with a 2020 first-pass
  analysis that used a different batch-correction approach.
- A CellChat-package-based cell-cell-communication exploration
  (`3.Result/CellChat*`, ~250 auto-generated per-pathway plots in the original
  project folder) is **not** included here; the methods draft
  (`manuscript/method.docx`) describes CellTalker as the algorithm used for
  the published receptor-ligand analysis, so CellChat appears to have been an
  exploratory alternative rather than what generated Figure 6c-d. If this call
  turns out to be wrong, the original files are still on the lab OneDrive.
- Two large, recently-dated summary tables (`monoMacroID.csv`,
  `overallCellIndentity.csv`, final cell-identity/cluster annotation) are kept
  alongside the source data on OneDrive rather than duplicated here.
