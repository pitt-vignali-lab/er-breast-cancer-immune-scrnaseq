# ER+ Breast Cancer Immune Landscape — scRNA-seq Analysis

Analysis code and figure-level results supporting the single-cell RNA-seq
sections of:

> Onkar S, Cui J, Zou J, Cardello C, Cillo AR, Uddin MR, Sagan A, Joy M,
> Osmanbeyoglu HU, Pogue-Geile K, McAuliffe PF, Lucas PC, Tseng GC, Lee AV,
> Bruno TC, Oesterreich S, Vignali DAA.
> **Immune landscape in estrogen receptor positive breast cancer reveals a
> divergent macrophage-driven microenvironment.**
> *Nature Cancer*. 2023;4(4):516–534.
> DOI: [10.1038/s43018-023-00527-w](https://doi.org/10.1038/s43018-023-00527-w) ·
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

This repository contains only the parts of the original analysis project that
fed into the **published main figures (Figure 4, Figure 5, Figure 6)**:
clustering/annotation of the CD45+ immune compartment, monocyte/macrophage
sub-clustering, RNA velocity/pseudotime, differential expression, SPIA pathway
enrichment, and macrophage–T cell receptor–ligand (cell–cell communication)
analysis. Flow cytometry, mIHC, and clinical/survival analyses (Figures 1–3, 7
and most Extended Data figures) were done outside this codebase and are not
included here.

Raw sequencing data (10x gene expression matrices, CITE-seq/hashing counts,
velocyto loom files) and large intermediate Seurat/CellChat `.rds`/`.RData`
objects (0.1–4 GB each) are **not** included — GitHub is not an appropriate
place for them. They remain on the lab's OneDrive project folder
(`Sayali/SayaliNew`). Scripts assume those objects exist locally in an
`1.Data/` / intermediate-object layout matching the original project.

## Repository layout

```
scripts/
  01_qc_clustering_annotation/     QC, initial clustering, SingleR/manual cluster
                                    labeling, Bhattacharyya distance, marker dotplot
                                    -> Figure 4a-d
  02_monocyte_macrophage_subclustering/
                                    Re-clustering of mono/macro cells, marker genes,
                                    DE, selected-gene heatmaps -> Figure 4e-f, 5a, 5d-e
  03_rna_velocity_pseudotime/       scVelo RNA velocity + latent time -> Figure 5b-c
  04_pathway_enrichment_SPIA/       SPIA pathway enrichment on cluster DE genes
                                    -> Figure 6a-b
  05_cell_cell_communication/       celltalker-based receptor-ligand analysis
                                    between monocytes/macrophages and T cells
                                    -> Figure 6c-d

results/
  figure4_immune_landscape_monocyte_macrophage/
  figure5_monomacro_trajectory_DE/
  figure6_pathway_cell_communication/
    each contains:
      source_data/   official per-panel source data (as deposited with the paper)
      figures/        rendered plots (png/pdf) produced by the scripts above
      *_tables/       supporting DE/GO/KEGG/receptor-ligand tables

environment/
  Analysis.Rproj                   RStudio project file used for the R scripts

draft/
  methods_draft.docx               internal methods draft used while writing the paper
```

## Notes on scope

- Scripts are the versions under the project's `2.Analysis/` pipeline, which
  were iterated on from 2021–2023 and are the versions that produced the
  `source_data/*.csv` files matching the published figure panels. An earlier
  2021 draft pass of the same pipeline (root-level `1.test.R`, `2.clustering.R`,
  `4.x`, `5.x` scripts) was superseded by this pipeline and is intentionally
  excluded.
- A `CellChat`-package-based cell-cell-communication exploration
  (`3.Result/CellChat*`, ~250 auto-generated per-pathway plots) also exists in
  the original project folder but does not correspond to the specific
  receptor-ligand circos plots shown in Figure 6c-d; it was left out as an
  exploratory alternative. Likewise a slingshot/diffusion-map trajectory
  exploration was superseded by the scVelo-based analysis actually shown in
  Figure 5b-c, and a GSEA pass (not cited in the Methods) was left out. If any
  of these calls are wrong, the original files are still on OneDrive
  (`Sayali/SayaliNew`) and can be added back.
- File names were normalized (figure-panel prefixes, English names, no spaces)
  for readability; original filenames are noted in code comments where they
  differ substantially.
