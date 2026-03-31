# Inverse Signal Importance (ISI)

A computational framework for estimating **time-varying signal importance** — how the influence of environmental signals on biological responses changes over time.

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
├── ISI.py                               # Core ISI module (EMAlgorithm class)
├── iaaft.py                             # IAAFT surrogate data generation
├── InverseSignalImportance_share.ipynb  # Main analysis notebook (ISI estimation)
├── RandomizeSignalImportance.ipynb      # Randomization test notebook (IAAFT surrogates)
├── Medaka_RNAseq_analysis_share.Rmd    # RNA-seq transcriptome analysis (R)
├── Medaka_RNAseq_analysis_share.html   # Knitted HTML output of the Rmd analysis
├── medaka_female_sex_hormone_evidence.html  # Reference table of female sex hormone-related genes
├── 20151001_20171015_env_inGSI.csv     # Environmental signals and GSI (gonadosomatic index) data
├── 20151001_20171015_raw_GSIdata_ci.csv # GSI measurements
├── org.Olatipes.eg.db/                 # Custom R annotation package for medaka (Oryzias latipes)
├── RNAseq_nakayama_2015to17/           # Transcriptome data (medaka brain, 2015–2017)
│   ├── group1/                         # Biological replicate 1
│   ├── group2/                         # Biological replicate 2
│   └── README.md
└── TV0603_alpha0win3noisemean/         # ISI model output
    └── SignalImportance_average.csv    # Estimated signal importance used in the downstream analysis
```

---

## Data

### Environmental and Reproductive Data

Outdoor measurements of medaka (*Oryzias latipes*) sampled over two years (October 2015 – October 2017).

| Column | Description |
|--------|-------------|
| `SR`   | Solar radiation |
| `WT`   | Water temperature |
| `DL`   | Day length (photoperiod) |
| `GSI`  | Gonadosomatic index (reproductive state proxy) |

### Transcriptome Data

Genome-wide transcriptome analysis of medaka brain regions (ventral telencephalon, hypothalamus, and pituitary) sampled under outdoor natural conditions.

- **GEO accession**: [GSE234401](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234401)
- Group 1 and Group 2 are biological replicates
- Low-quality sequencing reads were filtered prior to analysis.

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
    n_signals=3        # SR, WT, DL
)

rss_test, predictions = em.test(
    X_list=[X3],
    Y_list=[Y3],
    x_smooth=x_smooth,
    coef=coef
)
```

**Key parameters:**

| Parameter | Description |
|-----------|-------------|
| `window_size` | Number of consecutive days in each sliding input window (response function length) |
| `process_noise_var` | Controls how rapidly signal importance can change over time |
| `observe_noise_var` | Diagonal variance of observation noise (set from dGSI variance) |
| `observe_noise_cov` | Off-diagonal covariance of observation noise (set from dGSI covariance between replicates) |
| `alpha` | Ridge regularization strength |

### Surrogate Testing (`iaaft.py`, `RandomizeSignalImportance.ipynb`)

IAAFT surrogates are used to assess statistical significance of similarity between signal importance and temporal transcriptome profile. Each surrogate preserves the amplitude distribution and power spectrum of the original signal while destroying temporal correlations.

### Transcriptome Analysis (`Medaka_RNAseq_analysis_share.Rmd`)

Compares outdoor RNA-seq profiles to ISI-estimated signal importance trajectories. Identifies genes whose expression correlates with the time-varying importance of specific environmental signals.

- **`Medaka_RNAseq_analysis_share.html`**: Knitted HTML output for browsing the full analysis results without running R.
- **`medaka_female_sex_hormone_evidence.html`**: Reference table listing female sex hormone-related genes in medaka with roles in steroidogenesis/signaling, PubMed links, and KEGG Orthology entries. Used to annotate genes of interest in the transcriptome analysis.
- **`org.Olatipes.eg.db/`**: A custom R annotation package (OrgDb) for medaka (*Oryzias latipes*), providing gene-level mappings based on Entrez Gene identifiers. Required by `Medaka_RNAseq_analysis_share.Rmd` for gene annotation and GO/pathway analysis.

---

## Dependencies

### Python (version 3.12.11)

- `numpy`
- `pandas`
- `scipy`
- `scikit-learn`
- `tqdm`

Install with:

```bash
pip install numpy pandas scipy scikit-learn tqdm
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

`iaaft.py` is licensed under the **GNU Affero General Public License v3.0** (AGPL-3.0), following the original implementation by Bedartha Goswami.
https://github.com/mlcs/iaaft

Other files in this repository are provided for research reproducibility.
