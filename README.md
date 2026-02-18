[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18672348.svg)](https://doi.org/10.5281/zenodo.18672348)

# TFMethyl DNAm aging-clock

**Source code for the work "Enhancing the performance and interpretability of epigenetic clocks"**
🔗 [Preprint](https://www.biorxiv.org/content/10.1101/2025.10.07.680024v2)

---

## Overview

This repository contains the workflow, code, and model developed during the work. The project integrates both **R** and **Bash** code and provides **.Rmd** markdown for reproducible analysis.

---

## Repository Structure

| Element                  | Description                                                                                                                                                                                    |
| -----------------------  | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `./Run_TFMethyl`         | Standalone package implementation of the epigenetic age prediction model for humans. Contains the trained model, execution scripts, and a small toy dataset in the required input format.      |
| `./helper_scripts`       | Utility scripts for generating genome-wide TFBS annotations, performing footprinting analysis, and checking required R packages for `Analysis.Rmd`. See README inside for details.             |
| `Analysis.Rmd`           | R Markdown source used to generate figures and perform analyses for the study *“Enhancing the performance and interpretability of epigenetic clocks”*.                                         |
| `Analysis.html`          | Knitted HTML output of `Analysis.Rmd`, provided for convenient viewing in a web browser.                                                                                                       |

---

## 🧬 TFMethyl Usage

To run the TFMethyl clock, place a beta matrix file (named "beta_test_matrix.csv") at: `./Run_TFMethyl/` and run `source("./Run_TFMethyl/TFMethyl.R")` in an R environment. Detailed specifications for the input matrix, model, and usage are available at `./Run_TFMethyl/README`.

## 🔧 Analysis Reproducibility

`Analysis.Rmd` is dependent on bunch of input data structures mounted at `./Data/` directory, and installed CRAN/Bioconductor packages. Before running the script, do as follows:

```bash
wget https://zenodo.org/records/18672348/files/data.gz
tar -xzvf data.gz
```
Then first run `./helper_scripts/packages_check.R`, finally followed by the `Analysis.Rmd` in a R enviornment. 

## Notes

* The processed data (as`./Data`) required to run `Analysis.Rmd` - is the ONLY missing directory (~16GBs) in this repository (download from Zenodo). 
* All TFMethyl clock required data objects are included as (`./Run_TFMethyl/*.rds`). 

---
