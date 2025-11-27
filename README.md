# NonSyn Variant Inspector
### Shiny app for prioritizing non-synonymous variants using dbSNP + dbNSFP annotations

## Overview
NonSyn Variant Inspector is a Shiny application designed to prioritize non-synonymous variants by combining:
- GWAS variant tables (CHR, BP, P, SNP…)
- Automatic extraction of variants from genomic intervals via **bcftools**
- Full functional annotation using **dbNSFP 5.0a**
- Interactive visualizations: **Manhattan**, **Heatmaps**, **GO**, and **KEGG**
- Enriched tables with interactive links (ClinVar, dbSNP, GeneCards, UCSC)

The workflow consists of **4 preprocessing steps** plus a final visualization module.

---

## Project Structure
```
NonSyn-Variant-Inspector/
│
├── app.R                 
├── R/
│   ├── format_dbnsfp.R
│   └── clinical_metrics.R
├── README.md
└── resources/
```

---

# Requirements

## R Packages
Main dependencies:
- shiny, shinyjs, plotly, DT  
- dplyr, tidyr, readr  
- ggplot2, ggplotify, pheatmap  
- clusterProfiler, org.Hs.eg.db, GO.db  
- processx, conflicted, stringr  

---

# External Tools & Required Resources

## 1. **bcftools (required)**  
Used to extract variant sites from genomic intervals.

### Installation
**Conda (recommended):**
```bash
conda install -c bioconda bcftools
```

**Homebrew (macOS):**
```bash
brew install bcftools
```

**Official releases (source/binaries):**  
https://github.com/samtools/bcftools/releases

The app allows specifying the full path to the binary if it's not in `$PATH`.

---

## 2. **dbNSFP 5.0a (required)**  
Massive functional annotation dataset for non-synonymous variants.

### Download
**Official site:**  
dbNSFP is free for academic and non-commercial use under the CC BY-NC-ND 4.0 license.
Register with dbNSFP using your academic institutional email to receive the dbNSFP Access Code
Then use the Access Code to request a download for the current or past versions of dbNSFP. The download links will be sent to your registered email automatically.
[https://dbnsfp.houstonbioinformatics.org](https://www.dbnsfp.org)

You must extract and include in a single folder:
```
dbNSFP5.0a_variant.chr*.gz
dbNSFP5.0_gene.gz
search_dbNSFP50a.class
```

This folder is specified in **Step 3 — dbNSFP annotation**.

---

## 3. **dbSNP hg38 VCF (Broad Institute version, required)**  
Used in *Step 2* for extracting variants located in selected intervals.

### Required file:
```
resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz
```

### Download
Google Cloud (Broad Institute):

```
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi
```

This file is fully compatible with:
- bcftools
- hg38 genome
- UCSC Genome Browser
- Manhattan plot genome layout

---

# Workflow

## **Step 1 — Import GWAS Table**
Load a TSV/CSV containing at least:
- **CHR**, **BP**, **P**
- Optional: **SNP**, alleles, gene name, etc.

Generates:
- Manhattan plot  
- Table of significant SNPs  
- Auto-generated genomic intervals  

---

## **Step 2 — Extract Variants with bcftools**
Uses dbSNP VCF to extract:
- Only variants in selected genomic intervals  
- Auto-handling of `chr` prefix  
- Merging overlapping intervals  
- Optional VCFv4.0 forced header  

Output:
```
selected_sites.vcf
```

---

## **Step 3 — Annotate with dbNSFP**
Runs the Java dbNSFP search engine:

```
java -Xmx6g search_dbNSFP50a -i selected_intervals.vcf -o dbnsfp.out
```

---

## **Step 4 — Normalize dbNSFP Output**
Produces clean and structured files:
- `dbnsfp_normalized.csv`
- `dbnsfp_normalized.rds`

Cleans numeric fields, fixes separators, tags gene symbols, and standardizes chromosome fields.

---

# Visualization Modules

## **Heatmap View**
- Score / Rankscore toggle  
- Column scaling (0–1)  
- Chromosome filtering  
- Top-N SNP selection  
- PNG export  

## **dbNSFP Table**
- Functional class filtering (Pathogenicity, Conservation, LoF…)  
- GeneCards links  
- dbSNP links  
- ClinVar and OMIM fields  
- UCSC region link per SNP  

## **Manhattan & Genes**
- Multi-metric Manhattan plot  
- Supports all Score/Rankscore fields  
- Click → UCSC Genome Browser (centered + highlighted base)  
- Selectable UCSC mirror  

## **GO Enrichment**
- Biological Process, Cellular Component, Molecular Function  
- clusterProfiler + org.Hs.eg.db  
- Simulated FDR  
- PNG export  

## **KEGG Pathways**
- KEGG enrichment using `enrichKEGG`  
- FDR slider  
- Clickable KEGG pathway links  
- PNG export  

---

# Running the App

### Local:
```r
shiny::runApp("NonSyn-Variant-Inspector")
```

### RStudio:
Open `app.R` and click **Run App**.

---

# Required Input Files
- GWAS table (CSV/TSV)  
- dbSNP hg38 VCF (Broad Institute)  
- Full dbNSFP 5.0a folder  
- Java (for dbNSFP search engine)

---

# Notes
- Temporary directories are automatically deleted at the end of each Shiny session.
- Compatible with **macOS (Intel + Apple Silicon), Linux, and server installations**.
- dbNSFP lookup may take several minutes depending on system resources.

---

# Author
**Joan Fibla (jfibla)**  
Maintainer and developer.  
