## Overview

This repository contains the key scripts and custom R functions developed to support the analysis published in:

> **Adekoya AE, et al.** *Revealing Community Dynamics in Polymicrobial Infections through a Quantitative Framework.* ISME Communications, 2026.

The study uses a **quantitative accuracy scoring framework** to evaluate how closely in vitro polymicrobial infection models mimic in situ microbial physiology observed in human chronic wound infections 

---

## Background

Chronic wound infections are polymicrobial: multiple bacterial species coexist and interact in ways that drive disease persistence and severity. Building laboratory models that mimic these communities is critical for studying infection biology and testing interventions. However, the previously used quantitative framework for assessing how accurately models reflect actual human infection transcriptomes was based on a single-infection model.

This framework addresses that gap by extending the **Accuracy Scores (AS2)**  beyond a single-species infection model to a polymicrobial model. The AS2 framework uses a z-score-based metric to quantify, at the level of individual functional gene categories (COGs), how closely a model community's transcriptional profile matches a target human infection dataset. The framework includes leave-one-out cross-validation for self-validation and sunburst visualizations for intuitive representation of community-level accuracy.

---

## Repository Contents


accuracy_functions.R | Core custom R functions for AS2 scoring, z-score computation, leave-one-out validation, and sunburst visualization |
accuracy_score.sh | SLURM-based HPC submission script for running accuracy scoring at scale |
metaphlan4.sh | Pipeline script for taxonomic profiling using MetaPhlAn4 |
humann4.sh | Pipeline script for functional profiling using HUMAnN3/4 |

---

## Core Functions (accuracy_functions.R)

### score_target_vs_model_SA(counts_DF, Target_samplenames, Model_samplenames)
Computes per-gene z-scores comparing a model community's transcriptional profile to a target (human infection) dataset. Returns a data frame of mean z-scores (penalty scores) per COG functional category.

- **Input:** Normalized counts data frame; lists of target and model sample names
- **Output:** Data frame of COG-level penalty scores
- **Method:** Per-gene z-score normalization against target mean ± SD; mean penalty computed across model replicates

### score_target_vs_model(Target_samplenames, Model_samplenames)
Streamlined version operating on a globally defined normalized counts object. Designed for iterative scoring across multiple model comparisons.

### model_selfvalidation_leaveout(sample_DF1, left_out_number)
Leave-one-out cross-validation function. Randomly holds out a specified number of human infection samples as a test set and scores remaining samples as the training reference — enabling internal validation of the framework's consistency.

### draw_sunburst_STANDARD_SA(SA_df_counts, COG_SA, sb, target_list, model_list, maxlim, node_labs)
Generates a sunburst visualization of AS2 scores organized by COG functional hierarchy (Meta → Main categories). Scores are mapped to a YlGnBu color scale representing the proportion of genes within each category that fall within the accuracy threshold.

### draw_sunburst_NO_TEXT(...)
Identical to above with node labels suppressed — used for cleaner figure production in manuscript panels.

---

## Dependencies

```r
# R packages required
library(dplyr)
library(tidyverse)
library(cowplot)
library(zeallot)
library(devtools)
library(reshape2)
library(ggsunburst)   # devtools::install_github("didacs/ggsunburst")
library(readxl)
library(scico)
library(stringr)
library(purrr)
library(tibble)
library(RColorBrewer)
```

---

## Usage

```r
# Source the functions
source("accuracy_functions.R")

# Score a model against a target human infection dataset
results <- score_target_vs_model_SA(
  counts_DF = your_normalized_counts,
  Target_samplenames = c("human_sample_1", "human_sample_2", "human_sample_3"),
  Model_samplenames = c("model_sample_1", "model_sample_2")
)

# Run leave-one-out cross-validation
validation <- model_selfvalidation_leaveout(sample_DF1 = human_samples, left_out_number = 2)

# Generate sunburst visualization
sunburst_output <- draw_sunburst_STANDARD_SA(
  SA_df_counts = your_counts,
  COG_SA = cog_annotations,
  sb = sunburst_data,
  target_list = human_list,
  model_list = model_list,
  maxlim = 2,
  node_labs = TRUE
)
```

---

## Computational Environment

All analyses were run on R studio and Linux-based HPC clusters using SLURM job scheduling. See `accuracy_score.sh` for an example SLURM submission script.

---

## Citation

If you use these scripts in your work, please cite:

> Adekoya AE, et al. *Revealing Community Dynamics in Polymicrobial Infections through a Quantitative Framework.* ISME Communications, 2026.

---

## Contact

**Aanuoluwa Enitan Adekoya**
PhD Candidate, Microbiology — University of Tennessee, Knoxville
aadekoya@vols.utk.edu | [LinkedIn](https://linkedin.com/in/aanuoluwaadekoya) | [GitHub](https://github.com/Aanuoluwaduro)
