# seasonal-mangrove-bacterial-community
Analysis code for the study of seasonal variation in mangrove sediment bacterial communities.
# Seasonal Mangrove Bacterial Community Analysis

This repository contains the R scripts used for the statistical and ecological analyses in our study of seasonal variation in mangrove sediment bacterial communities.

## Scripts

### `alpha_environment_PERMANOVA.R`

This script performs analyses of bacterial alpha diversity, environmental variables, and community composition.

Main analyses include:

* Linear models for bacterial richness, Chao1, Shannon diversity, and Faith's phylogenetic diversity
* Effects of pH after accounting for season, site, sediment compartment, and other measured environmental variables
* Variance inflation factor (VIF) analysis
* Bray-Curtis dissimilarity calculation
* PERMANOVA
* Partial constrained ordination
* Multivariate dispersion analysis
* Sensitivity analysis using residual pH variation

Main input files:

* `envandalpha1.txt`
* `ASV_table.txt`

### `Partial_CAP.R`

This script performs principal coordinates analysis (PCoA) and partial constrained analysis of principal coordinates (CAP) based on Bray-Curtis dissimilarity.

Main analyses include:

* Hellinger transformation
* Bray-Curtis dissimilarity
* PCoA with Lingoes correction
* Species projection onto ordination space
* PCoA visualization by season and sediment compartment
* Partial CAP analysis while accounting for site and sediment compartment

Main input files:

* `asv_table.txt`
* `group.txt`

### `SEM_piecewise.R`

This script performs piecewise structural equation modeling to evaluate relationships among sediment pH, bacterial richness, homogeneous selection, and the richness-difference component of beta diversity.

Main analyses include:

* Conversion of pairwise ecological-process metrics to sample-level variables
* Piecewise structural equation modeling
* Model fit assessment using Fisher's C
* Standardized path coefficients
* R-squared calculation
* Variance inflation factor analysis
* Breusch-Pagan tests
* Influence diagnostics
* Sensitivity analyses using median aggregation and within-site sample pairs

Main input files:

* `icamp.txt`
* `rep.txt`
* `asv_table.txt`
* `env.txt`

### `MPD_panel_f.R`

This script evaluates the phylogenetic clustering of selected phylogenetic bins using mean pairwise phylogenetic distance (MPD).

Main analyses include:

* Matching phylogenetic bins to the phylogenetic tree
* Calculation of observed MPD
* Generation of a null MPD distribution using 9,999 randomizations
* One-sided and two-sided permutation tests
* Standardized effect size calculation
* Visualization of the observed MPD relative to the null distribution

Main input files:

* `datanet1.txt`
* `root.tree`

## Software

The analyses were performed in R.

Main R packages used in these scripts include:

* `vegan`
* `adespatial`
* `ape`
* `data.table`
* `tidyverse`
* `dplyr`
* `tidyr`
* `ggplot2`
* `piecewiseSEM`
* `car`
* `lmtest`
* `permute`

Package versions and additional session information can be obtained using `sessionInfo()` in R.

## Usage

Download or clone this repository and place the required input files in the corresponding working directory.

Before running each script, modify the working directory specified by:

```r
setwd("your_path")
```

Then run the corresponding R script.

For example:

```r
source("SEM_piecewise.R")
```

Input file names and required column names are specified within each script.

## Data availability

The sequencing reads generated in this study have been deposited in the NCBI database under BioProject accession number **PRJNA1420879**.

The analysis scripts used in the study are publicly available in this repository.

## License

This repository is distributed under the terms specified in the `LICENSE` file.
