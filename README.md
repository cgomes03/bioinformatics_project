# Application of rMTA-EA Pipeline to Multi-Gene Metabolic Optimization in Neurodegenerative Diseases

This repository contains the implementation of a robust computational pipeline that combines the Robust Metabolic Transformation Algorithm (rMTA) with Evolutionary Algorithms (EAs) to identify effective multi-gene knockout strategies in the context of neurodegenerative diseases, with a focus on Alzheimer’s Disease (AD).

## Overview

The pipeline integrates transcriptomic data into genome-scale metabolic models (GEMs), in specific Recon3D, reconstructs context-specific metabolic networks, and applies rMTA to score candidate interventions based on the robust transformation score (rTS). Evolutionary Algorithms are then used to optimize combinations of gene knockouts .

## Structure

```
BIOINFORMATICS_PROJECT/
├── Article and Presentation/
│   ├── catarina_gomes_pg55694_pt1.pdf
│   └── presentation_pg55694.pdf
├── data/
│   ├── Context Specific/              # iMAT-generated tissue-specific models
│   ├── GSE203206/                     # GEO transcriptomic data
│   └── Recon3D/                       # Recon3D, GEM and derived analysis data and Samples ACHR
├── results/
│   ├── outputs/                       # Outputs of individual EA+rMTA runs
│   └── result_merged/                 # Merged analysis, visualizations, pruning
├── rMTA_pipeline/
│   ├── analysis/                      # Jupyter notebooks for result interpretation
│   ├── dge_r/                         # R scripts for DGE processing
│   ├── dict/, differential_expression/, environments/, expression_data/
│   ├── scripts/                       # Data conversion and sampling
│   └── utils/                         # Helper scripts and outputs
├── src/
│   └── EA_rMTA/                       # Core code for EAs, rMTA and Simplification 
└── README.md
```

## Dataset

Transcriptomic data from [GEO: GSE203206](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203206), which profiles Alzheimer's patients (early- and late-onset) and healthy controls. Data is integrated into the Recon3D GEM to generate condition-specific models.

## Publications and Acknowledgements

This work was developed as part of the curricular unit 'Project in Bioinformatics' by Catarina Gomes at the University of Minho, under the supervision of Bruno Sá and Miguel Rocha.