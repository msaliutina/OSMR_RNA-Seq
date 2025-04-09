# OSMR RNA-Seq Analysis

This repository contains the code used for RNA-Seq differential expression analysis of **OSMR** (Oncostatin M Receptor) expression across different sample groups, as presented in the Nature Immunology publication (2025):  
**The interleukin 22–oncostatin M axis promotes intestinal inflammation and tumorigenesis**

## Overview

The analysis includes the following steps:

- Preprocessing of bulk and single cell RNA-Seq data
- Exploratory data analysis
- Generation of visualizations (PCA, heatmaps, MA plots, volcano plots)
- cell-cell interaction analysis
- etc
  
## Repository Structure
```
OSMR_RNA-Seq/
├── data/ # Input data and intermediate objects
│ ├── bulk_data/ # Bulk RNA-Seq count matrices and metadata
│ ├── cell_lines/ # Processed data for cell line analysis
│ ├── cellphonedb_output/ # Output from CellPhoneDB interaction analysis and futher visualization
│ ├── corr_analysis/ # Data for correlation analysis of UNIFI/IBDome/Mount Sinai cohorts
│ ├── sc_final/ # Final single-cell Seurat objects (download required, check "Data availability")
│ └── sc_raw/ # Seurat objects before preprocessing (download required, check "Data availability")
├── preprocessing/ # Input data and intermediate objects
│ ├── preprocessing_EpCam.R # EpCam analysis
│ ├── preprocessing_LPMC.R # LPMC analysis
│ ├── preprocessing_Stroma.R # Stroma analysis
├── Ext_1.R                   # Extended Data Figure 1
├── Ext_2_CellPhoneDB.R       # CellPhoneDB ligand-receptor interaction analysis
├── Ext_6_ibdome.R            # Extended Data Figure 6 — IBDome cohort
├── Ext_6_sinai.R             # Extended Data Figure 6 — Sinai cohort
├── Ext_6_UNIFI.R             # Extended Data Figure 6 — UNIFI cohort
├── Ext_9.R                   # Extended Data Figure 9
├── Fig1.R                    # Main Figure 1
├── Fig2.R                    # Main Figure 2
├── Fig6_7.R                  # Main Figures 6 and 7
├── FigS1.R                   # Supplementary Figure S1
├── FigS2.R                   # Supplementary Figure S2
├── FigS3.R                   # Supplementary Figure S3
├── FigS4.R                   # Supplementary Figure S4
└── README.md                 # Project documentation
```

## Data Availability

Seurat objects which were taken for the further preprocessing are available at the following link:

[Download "raw" Seurat objects](https://drive.google.com/drive/folders/1uh_MJCUf4giD7afw03oiWz9x-1dS5Tjr?usp=share_link)

Finalized versions of seurat objects are available here:

[Finalized Seurat objects](https://drive.google.com/drive/folders/1NZb-0ivlP_sh3VaOPZODYFo1bdnwDlQz?usp=share_link)




