# Systematic functional screening of chromatin factors identifies strong lineage and disease dependencies in normal and malignant haematopoiesis.
This folder contains the R code to analyze FACS-based and Perturb-seq CRISPR screens. R scripts in this folder are organized by the type of analysis.

To run this code you need the Singularity container, the data processed by cell ranger, and then use renv to install the packages in the file renv.lock.

## Contents
- Scripts to obtain external data, start with EXT
- Scripts to analyze FACS-based data start with POOLED
- Scripts to analyze Perturb-seq data start with SCRNA
    - Basic analysis (QC, integration, duplet detection, cell type identification, analysis of marker genes) in files SCRNA_01 to SCRNA_07
    - Projection of datasets in files SCRNA_08
    - Results of the above analyses are collected in SCRNA_10
    - Cluster and cell type analyses of the above are in files SCRNA_20 to SCRNA_22
    - Differential expression analyses in files SCRNA_30 to SCRNA_41
- Additional functions are defined in files starting with FUNC
- Publication figures are produced in files starting with FIG
    - Analyses of FACS-based screens in file FIG_01
    - Analyses of Perturb-seq data focused on UMAPs, clusters, and cell types in FIG_02
    - Analyses of Perturb-seq data focused on differential expression analyses in FIG_03
    - Exporting of supplementary tables in FIG_80
