# rgm: Random Graphical Models for Network Analysis

`rgm` is an R package that implements state-of-the-art Random Graphical Models (RGM) for the analysis of complex network data. It excels in handling heterogeneous data across various environments, offering a powerful tool for exploring intricate network interactions and structural relationships.

## Key Features

- **Joint Inference Across Multiple Environments**: `rgm` enables simultaneous analysis of data from diverse environments or communities, providing a comprehensive understanding of complex network interactions.
- **Advanced Random Graphical Modeling**: The package employs sophisticated RGM methodologies to effectively handle heterogeneity and quantify structural relationships between datasets.
- **Integration of External Covariates**: Users can incorporate external covariates at both node and interaction levels, allowing for a more nuanced analysis of network data.
- **Bayesian Framework**: `rgm` uses a Bayesian approach to robustly quantify parameter uncertainty, offering more reliable analytical results.
- **Network Posterior Visualization and Analysis**: The package includes tools for visualizing and interpreting network posteriors, revealing core structures and elucidating individual differences and relationships.

## Installation

Install the latest version of `rgm` from GitHub using the following commands in R:

```R
install.packages("devtools")
devtools::install_github("franciscorichter/rgm")
```

## Usage

For detailed instructions on using `rgm` for data analysis, refer to the package vignette and documentation:

```R
library(rgm)
vignette("rgm")
```

**Note**: While initially designed for microbiome analysis, `rgm` is broadly applicable across various fields requiring advanced graphical modeling of network data.

## Principal Reference

The methodologies implemented in the rgm package are principally derived from the work described in Vinciotti, V., Wit, E., & Richter, F. (2023). "Random Graphical Model of Microbiome Interactions in Related Environments." arXiv preprint arXiv:2304.01956.


