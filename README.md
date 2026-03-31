# Inverse Signal Importance (ISI)

A computational framework for estimating **time-varying signal importance**: how the influence of environmental signals on biological responses changes over time.

This repository provides the analysis code for identifying which environmental cues (photoperiod, water temperature, solar radiation) drive seasonal reproductive changes in the medaka (*Oryzias latipes*) under natural outdoor conditions.

---

## Overview

Environmental signals often vary in their importance to biological systems across seasons. Standard regression approaches assume fixed coefficients, but the relative importance of signals (e.g., day length vs. temperature) may shift throughout the year.

**ISI** estimates this time-varying importance by combining:

- A **continuous-time Kalman filter** and **RTS smoother** to track latent signal importance states
- **Ridge regression** to estimate response functions mapping environmental signals to biological outputs
- A **hybrid EM algorithm** that iterates between these two steps until convergence

---

## Repository Structure

```
.
+-- ISI.py                              # Core ISI module (EMAlgorithm class)
+-- InverseSignalImportance.ipynb       # Main analysis notebook (ISI estimation)
+-- README.md
+-- data/
|   +-- 20151001_20171015_env_inGSI.csv      # Environmental signals and GSI data
|   +-- 20151001_20171015_raw_GSIdata_ci.csv # Raw GSI measurements
|   +-- RNAseq_nakayama_2015to17/
|       +-- group1/                          # Biological replicate 1
|       +-- group2/                          # Biological replicate 2
+-- iaaft/
|   +-- iaaft.py                             # IAAFT surrogate data generation
|   +-- RandomizeSignalImportance.ipynb      # Randomization test notebook
+-- results/
|   +-- SignalImportance_average.csv         # Estimated signal importance
+-- rnaseq_analysis/
    +-- Medaka_RNAseq_analysis.Rmd           # RNA-seq transcriptome analysis (R)
    +-- Medaka_RNAseq_analysis_share.html    # Knitted HTML output
    +-- medaka_female_sex_hormone_evidence.html
    +-- pnas2313514120sd07.csv               # Circannual gene list (Nakayama et al. 2023)
    +-- org.Olatipes.eg.db/                  # Custom R annotation package for medaka
```

---

## Data

### Environmental and Reproductive Data

Outdoor measurements of medaka (*Oryzias latipes*) sampled over two years (October 2015 - October 2017). Files are located in `data/`.

| Column | Description |
|--------|-------------|
| SR | Solar radiation |
| WT | Water temperature |
| DL | Day length (photoperiod) |
| GSI | Gonadosomatic index (reproductive state proxy) |

### Transcriptome Data

Genome-wide transcriptome analysis of medaka brain regions (ventral telencephalon, hypothalamus, and pituitary) sampled under outdoor natural conditions. Raw data are located in `data/RNAseq_nakayama_2015to17/`.

- GEO accession: [GSE234401](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234401)
- Group 1 and Group 2 are biological replicates
- Low-quality sequencing reads were filtered prior to analysis

---

## Methods

### ISI Module (`ISI.py`)

The `EMAlgorithm` class implements the core algorithm:

```python
import ISI

em = ISI.EMAlgorithm(
    window_size=3,
    process_noise_var=10.0,
    observe_noise_var=1.5,
    observe_noise_cov=0.35,
    alpha=0.0
)

x_smooth, P_smooth, coef, rss, rss_history = em.fit(
    X_list=[X1, X2],   # Sliding-window input matrices (one per replicate group)
    Y_list=[Y1, Y2],   # Observed dGSI time series
    times=time_array,
    n_signals=4        # SR, WT, DL, GSI
)

rss_test, predictions = em.test(
    X_list=[X3],
    Y_list=[Y3],
    x_smooth=x_smooth,
    coef=coef
)
```

Key parameters:

| Parameter | Description |
|-----------|-------------|
| `window_size` | Number of consecutive days in each sliding input window |
| `process_noise_var` | Controls how rapidly signal importance can change over time |
| `observe_noise_var` | Diagonal variance of observation noise (set from dGSI variance) |
| `observe_noise_cov` | Off-diagonal covariance of observation noise (set from dGSI covariance between replicates) |
| `alpha` | Ridge regularization strength |

### Main Analysis Notebook (`InverseSignalImportance.ipynb`)

Runs the full ISI pipeline in the following order:

1. **Data Preprocessing**: loading and interpolation of GSI time series, standardization of variables, construction of sliding-window input matrix
2. **ISI: Single-Group Demonstration**: application to Data Group 1, visualization of RSS, response function, and signal importance
3. **Hyperparameter Selection via Grid Search**: grid search over `window_size` and `process_noise_var`
4. **Full Model: All Three Data Groups**: full pipeline across all groups, aggregated results as reported in the paper

Estimated signal importance is saved to `results/SignalImportance_average.csv`. The notebook also launches `iaaft/RandomizeSignalImportance.ipynb` as a subprocess for surrogate testing.

### Surrogate Testing (`iaaft/iaaft.py`, `iaaft/RandomizeSignalImportance.ipynb`)

IAAFT surrogates are used to assess statistical significance of similarity between signal importance and temporal transcriptome profiles. Each surrogate preserves the amplitude distribution and power spectrum of the original signal while destroying temporal correlations.

Surrogate results are saved to `results/` as `randomized_SRimp.csv`, `randomized_WTimp.csv`, `randomized_DLimp.csv`, and `randomized_GSIimp.csv`.

### Transcriptome Analysis (`rnaseq_analysis/Medaka_RNAseq_analysis.Rmd`)

Compares outdoor RNA-seq profiles to ISI-estimated signal importance trajectories. Identifies genes whose expression correlates with the time-varying importance of specific environmental signals.

- **`Medaka_RNAseq_analysis_share.html`**: Knitted HTML output for browsing results without running R
- **`medaka_female_sex_hormone_evidence.html`**: Reference table of female sex hormone-related genes in medaka (steroidogenesis/signaling, PubMed links, KEGG Orthology entries)
- **`pnas2313514120sd07.csv`**: Circannual gene list from Nakayama et al. 2023
- **`org.Olatipes.eg.db/`**: Custom R annotation package (OrgDb) for medaka (*Oryzias latipes*), required for gene annotation and GO/pathway analysis

The Rmd reads data from `../data/RNAseq_nakayama_2015to17/` and signal importance from `../results/SignalImportance_average.csv`, relative to the `rnaseq_analysis/` directory.

---

## Dependencies

### Python (version 3.12.11)

- `numpy`, `pandas`, `scipy`, `scikit-learn`, `tqdm`, `matplotlib`, `seaborn`

```
pip install numpy pandas scipy scikit-learn tqdm matplotlib seaborn
```

### R (version 4.4.1)

- `tidyverse`, `stringr`, `gplots`, `matrixStats`
- `extrafont`, `ggpubr`, `gt`, `readxl`, `curl`
- `ggVennDiagram`, `pbapply`, `reticulate`

---

## Citation

If you use this code or data, please cite the associated publication (details to be added upon publication).

---

## License

`iaaft/iaaft.py` is licensed under the [GNU Affero General Public License v3.0](https://www.gnu.org/licenses/agpl-3.0.html) (AGPL-3.0), following the original implementation by Bedartha Goswami (https://github.com/mlcs/iaaft).

Other files in this repository are provided for research reproducibility.
