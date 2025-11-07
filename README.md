# RNA-seq Differential Expression Analysis Pipeline

A comprehensive, reproducible pipeline for bulk RNA-seq analysis from raw FASTQ files to publication-ready results. This pipeline was developed for analyzing iPSC-derived neural crest cells in a Treacher Collins syndrome study.

[![Language](https://img.shields.io/badge/Language-R%20%7C%20Bash-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## 📋 Overview

This pipeline performs end-to-end RNA-seq analysis including:
- Quality control assessment (MultiQC)
- Read alignment and quantification (STAR + RSEM)
- Differential expression analysis (DESeq2)
- Publication-quality visualizations

**Study Design:** Patient vs. Control comparison in iPSC-derived cells
- **Controls:** 3 biological replicates (2 conditions: treated/untreated)
- **Patients:** 5 biological replicates (2 conditions: treated/untreated)
- **Total samples:** 16 paired-end RNA-seq samples

## 🔧 Features

- **Automated workflow** from FASTQ to results
- **Robust normalization** using DESeq2 with gene length correction
- **Quality control** at multiple steps
- **Statistical rigor** with FDR correction and effect size filtering
- **Reproducible** with detailed documentation
- **Publication-ready figures** (PCA, MA plots, heatmaps)

## 📊 Example Results

### Principal Component Analysis
![PCA Plot](results/figures/PCA_plot.png)

### MA Plot - Differential Expression
![MA Plot](results/figures/MA_plot.png)

### Expression Heatmap
![Heatmap](results/figures/expression_heatmap.png)

## 🛠️ Requirements

### Software Dependencies
- **MultiQC** v1.14+
- **STAR** v2.7.11b
- **RSEM** v1.3.3+
- **R** v4.2+ with packages:
  - DESeq2
  - tximport
  - ggplot2
  - pheatmap
  - tidyverse
  - See `environment/requirements.txt` for complete list

### Hardware Recommendations
- **RAM:** 32GB minimum (64GB recommended)
- **CPU:** 8+ cores
- **Storage:** ~50GB for human genome reference + data

## 📦 Installation

### Clone the repository
```bash
git clone https://github.com/yourusername/RNA-seq-differential-expression-pipeline.git
cd RNA-seq-differential-expression-pipeline
```

### Install R dependencies
```r
# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "tximport", "apeglm"))

# Install CRAN packages
install.packages(c("tidyverse", "pheatmap", "RColorBrewer", 
                   "ggrepel", "GGally"))
```

### Download genome reference (human GRCh38)
```bash
# GENCODE v45 annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.primary_assembly.annotation.gtf.gz
gunzip gencode.v45.primary_assembly.annotation.gtf.gz

# GENCODE v45 genome
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
```

## 🚀 Quick Start

### 1. Prepare your data
Place raw FASTQ files in `data/raw/` and create a sample metadata file:

```bash
# data/raw/samples.txt
run,Condition
Sample1_R1.fastq,Control
Sample2_R1.fastq,Patient
```

### 2. Run quality control
```bash
bash scripts/01_quality_control.sh
```

### 3. Align and quantify reads
```bash
# Generate STAR+RSEM reference (one-time setup)
bash scripts/02_alignment_quantification.sh --prepare-reference

# Run alignment and quantification
bash scripts/02_alignment_quantification.sh --run-samples
```

### 4. Preprocess for DESeq2
```bash
bash scripts/03_preprocessing.sh
```

### 5. Run differential expression analysis
```r
source("differential_expression.R")
```

## 📁 Pipeline Steps

### Step 1: Quality Control
- Runs FastQC on all samples
- Generates MultiQC summary report
- **Output:** `results/qc_reports/multiqc_report.html`

### Step 2: Alignment & Quantification
- Aligns reads to genome using STAR
- Quantifies gene expression with RSEM
- **Output:** `*.genes.results` files

### Step 3: Preprocessing
- Fixes zero-length genes (DESeq2 compatibility)
- Prepares data for import
- **Output:** `*_modificado.genes.results`

### Step 4: Differential Expression
- Imports data with tximport
- Normalizes counts (median-of-ratios)
- Performs statistical testing (Wald test)
- Applies log fold-change shrinkage (apeglm)
- Generates visualizations

**Key Parameters:**
- Significance threshold: adjusted p-value < 0.005
- Effect size: log2FC > 0 (up) or < 0 (down)
- Multiple testing correction: Benjamini-Hochberg FDR

## 📈 Output Files

```
results/
├── figures/
│   ├── PCA_plot.png
│   ├── correlation_heatmap.png
│   ├── dispersion_plot.png
│   ├── MA_plot.png
│   └── expression_heatmap.png
├── tables/
│   ├── res_table_ipsc_cont_pat.txt          # All results
│   ├── ipsc_up_pat_cont.txt                 # Upregulated genes
│   └── ipsc_down_pat_cont.txt               # Downregulated genes
└── qc_reports/
    └── multiqc_report.html
```

## 🔬 Analysis Details

### Normalization Strategy
- **Method:** DESeq2 median-of-ratios normalization
- **Gene length correction:** Handled by tximport (RSEM effective lengths)
- **Zero-length genes:** Replaced with length=1 (recommended by M. Love)

### Statistical Testing
- **Method:** DESeq2 Wald test
- **Shrinkage:** apeglm for effect size estimation
- **Multiple testing:** Benjamini-Hochberg FDR
- **Significance:** svalue (s-value from lfcShrink) < 0.005

### Quality Control Metrics
- Sample correlation matrix
- Principal component analysis
- Dispersion estimates
- MA plots

## 📊 Typical Results Summary

For this dataset:
- **Total genes tested:** ~60,000
- **Genes passing filters:** ~44,000
- **Significantly differentially expressed:** ~X genes (adjust based on your results)
  - Upregulated in patients: ~Y genes
  - Downregulated in patients: ~Z genes

## 🤝 Use Cases

This pipeline is suitable for:
- Bulk RNA-seq differential expression analysis
- Case-control study designs
- iPSC-derived cell comparisons
- Disease mechanism studies
- Pathway enrichment analysis (preprocessing)

## 📖 Citation

If you use this pipeline, please cite:

```
Avelino, Camila. (2024). RNA-seq Differential Expression Analysis Pipeline. 
GitHub repository: https://github.com/yourusername/RNA-seq-differential-expression-pipeline
```

**Key publications this pipeline is based on:**
- Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*.
- Dobin, A., et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*.
- Li, B., & Dewey, C.N. (2011). RSEM: accurate transcript quantification from RNA-Seq data. *BMC Bioinformatics*.

## 🐛 Troubleshooting

### Common Issues

**Problem:** "Error in tximport: gene lengths are zero"
```r
# Solution: Run preprocessing script to replace zero-length genes with 1
bash scripts/03_preprocessing.sh
```

**Problem:** STAR alignment fails with "not enough memory"
```bash
# Solution: Reduce the number of threads or increase the available RAM
--limitBAMsortRAM 30000000000  # Add to STAR command
```

**Problem:** DESeq2 fails with "matrix is not full rank"
```r
# Solution: Check your design matrix for confounding variables
# Ensure you have biological replicates for each condition
```

## 📧 Contact

**Camila Avelino**  
PhD, Genetics  

For questions, suggestions, or collaboration:
- GitHub: [@ccavelino](https://github.com/cavelino)
- Email: avelino.c@icloud.com
- LinkedIn:(https://www.linkedin.com/in/camila-c-avelino/))

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


---

**⭐ If you find this pipeline useful, please consider giving it a star!**

**🔧 Need help with your RNA-seq analysis?** I'm available for bioinformatics consulting and data analysis projects. Contact me for rates and availability.
