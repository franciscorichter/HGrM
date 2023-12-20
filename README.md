# rgm: Advanced Random Graphical Models for Network Analysis

`rgm` is an R package that implements advanced Random Graphical Models (RGM) for the analysis of complex network data. It is particularly adept at handling heterogeneous data across various environments, making it a versatile tool for exploring intricate network interactions and structural relationships.

## Features

### Joint Inference Across Multiple Environments
`rgm` facilitates the simultaneous analysis of data from diverse environments or communities, enabling a more comprehensive understanding of complex network interactions.

### Advanced Random Graphical Modeling
The package utilizes sophisticated RGM methodologies to manage heterogeneity and quantify structural relationships between different datasets.

### Incorporation of External Covariates
Users can integrate external covariates at both the node and interaction levels, allowing for a richer and more nuanced analysis of network data.

### Bayesian Framework
`rgm` employs a Bayesian approach, ensuring robust quantification of parameter uncertainty and providing more reliable analytical results.

### Network Posterior Visualization and Analysis
The package offers tools for visualizing and interpreting network posteriors, revealing stable core structures, and elucidating individual differences and relationships.

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
vignette("rgm-introduction")
```

**Note**: While initially designed for microbiome analysis, `rgm` is broadly applicable across various fields requiring advanced graphical modeling of network data.

## Citation

If you use `rgm` in your research, please cite:

Vinciotti, V., Wit, E., & Richter, F. (2023). Random graphical model of microbiome interactions in related environments. arXiv preprint arXiv:2304.01956.


