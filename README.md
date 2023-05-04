# Description 

rgm is an R package designed for the analysis of microbiota systems from metagenomic data across multiple body sites using Random Graphical Models (RGM). It captures the heterogeneity and structural relationships among different microbial communities and allows for the integration of external covariates at both the microbial and interaction levels, providing a comprehensive and adaptable approach to analyzing the complexity of microbiome data.

# Features:

Joint inference of microbiota systems: rgm enables the simultaneous analysis of metagenomic data from various body sites, allowing researchers to better understand the complex network of interactions among microbes.

Random Graphical Models (RGM): The package employs an RGM approach to account for heterogeneity across different environments while quantifying their structural relatedness.

Inclusion of external covariates: rgm allows users to incorporate external covariates at both the microbial and interaction levels, which further adapts to the richness and complexity of microbiome data.

Bayesian implementation: The Bayesian implementation of the RGM fully quantifies parameter uncertainty, providing robust and reliable results.

Network posterior analysis: rgm offers tools to visualize and interpret the microbiome network posteriors, revealing not only a stable core structure but also individual differences between various body sites and relationships between different classes of microbes.

    Taxonomical classification support: The package is capable of handling taxonomical classification data, allowing researchers to easily compare the structural similarity of microbial communities across body sites.

# Installation:

To install the latest version of rgm, use the following command in R:

install.packages("rgm")

# Usage:

For detailed instructions on how to use rgm for your metagenomic data analysis, refer to the package vignette and documentation:

```
library(rgm)
vignette("rgm-introduction")
```
