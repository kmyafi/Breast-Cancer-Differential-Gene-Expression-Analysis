# Breast Cancer Differential Gene Expression Analysis
Differential Gene Expression (DGE) Analysis in Microarray Data of Breast Cancer Subtypes

> **Objectives :**
> 1. Perform gene filtering to filter out less expressive genes.
> 2. Give recommendations for genes that are expressed differently between four cancer subtypes: Luminal A, Luminal B, Basal, and HER.
> 3. Investigate genes that are expressed differently between normal cell and cancer cell.

## Set
- **Dataset:** GSE45827 on breast cancer gene expression from [Curated Microarray Database (CuMiDa)](https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida/);
- **Platform:** GPL 570
- **Filtering:** `genefilter`;
- **Gene Expression Analysis:** `limma`, `t-test`;
- **Gene Ontology:** `GO.db`.

## About The Data
Dataset GSE45827 on breast cancer gene expression from CuMiDa
* 6 classes
* 54676 genes
* 151 samples

## Algorithm included
- [x] Sampling
- [x] Gene filtering
- [x] Gene expression analysis: *LIMMA (linear models for microarray data)*, *t*-test
- [x] Gene ontology
