# SGDP Diversity Analysis

This repository contains code, data inputs, and figure generation scripts for estimating and visualizing human genetic diversity across global ancestry categories using the Simons Genome Diversity Project (SGDP).

## 🔬 Overview

We compute nucleotide diversity (π) from SGDP chromosome 22 data, group populations into GWAS ancestry categories, and compare genetic diversity proportions with global census-based representation. Our analysis demonstrates that census size is a poor proxy for genomic diversity and offers data-driven insight into underrepresented populations.

## 📁 Repository Structure

- **DATA/**: SGDP variant and metadata files
- **PYTHON/**: Analysis and plotting scripts
- **FIGS/**: Final figures for submission (Figure 1a, 1b)
- **RESULTS/**: Text-based output summaries

## 🔧 Key Scripts

- `01-SGDP-pi-bootstrap-superpops.py`: Computes π estimates and bootstraps per ancestry group.
- `02-render-ancestry-proportions-latam.py`: Plots regional ancestry proportions across Latin America.

## 📊 Figures

- **Figure 1a**: Census population vs. genetic diversity by GWAS ancestry group
- **Figure 1b**: Ancestry composition of 19 Latin American countries (based on Adhikari et al., 2017)

## 📄 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## 📣 Citation

If you use this code, please cite:

Corpas M, et al. *Why Genomic Diversity Should Not Be Framed by Census Alone.* (Correspondence, under review)

