README
================

# Population Viability Analysis of Fishing Cats (*Prionailurus viverrinus*)

## Overview

This repository accompanies the manuscript:

> \[Insert manuscript citation\]

The repository contains metadata and R scripts necessary to reproduce
the demographic analyses, stochastic population simulations, and
graphical outputs presented in the manuscript.

The archived repository is available through Zenodo:

**DOI:** <https://doi.org/XXXXXX>

## Repository Contents

Four R script files are provided.

### 1. Population Viability Analysis:

This file contains three matrix population models developed to estimate
population viability.

- **Model C1:** Two-year population model (population growth rate only)
- **Model C2:** Two-year population model incorporating vital rates
- **Model C3:** One-year population model incorporating annual vital
  rates

### 2. Population Viability Analysis: Scenario-Based Simulations

This file contains simulation models used to evaluate population
viability under alternative demographic scenarios.

- **Model 1:** Age-structured vital rate model (baseline scenario)
- **Model 2:** Management scenario 1 (5% reduction in survival rates)
- **Model 3:** Management scenario 2 (5% reduction in recruitment rates)

### 3. Asymptotic Sensitivity Analysis: Matrix Elements

This file contains code for calculating:

- Matrix element sensitivities
- Matrix element elasticities

### 4. Asymptotic Sensitivity Analysis: Vital Rates

This file contains code for calculating:

- Vital rate sensitivities
- Vital rate elasticities

## Associated Data and Metadata

Datasets and associated metadata used in this study originated from a
long-term fishing cat monitoring programme and are described in detail
in the accompanying manuscript:

> \[Insert manuscript citation\]

Users are strongly encouraged to read the manuscript before using or
interpreting the scripts provided in this repository.

## Reproducibility

Analyses were conducted in R.

To reproduce the analyses:

1.  Download all repository files.
2.  Open the corresponding R script.
3.  Run scripts sequentially as described in the manuscript.
4.  Generated outputs should reproduce the reported population growth
    rates, sensitivity analyses, elasticity analyses, and simulation
    outcomes.

## Citation

Users of the data files are encouraged to cite both the data package and
the associated data paper. The citation of the data package is:

### Repository Citation

> \[Zenodo repository citation\]

### Manuscript Citation

> [Manuscript citation](#manuscript-citation)

## Contact and Collaboration

We encourage parties interested in re-using the codes to contact the
paper authors KP, MJG, CRN or DN for a discussion of research ideas and
potential opportunities for collaboration.

### Contact Authors

Contact details are available in the associated manuscript.
