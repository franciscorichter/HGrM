# Description:
rgm is an R package designed to provide a comprehensive and adaptable framework for analyzing complex network data using Random Graphical Models (RGM). This package can account for the heterogeneity and structural relationships among different communities or networks, and it allows users to incorporate external covariates at both the node and interaction levels.

# Features:

    Joint Inference of Systems: rgm enables the simultaneous analysis of data from various fields, allowing for a deeper understanding of complex network interactions.

    Random Graphical Models: The package employs RGM to handle heterogeneity across different environments and quantify structural relatedness.

    Inclusion of External Covariates: Users can incorporate external covariates at both the node and interaction levels to accommodate the richness and complexity of network data.

    Bayesian Implementation: This feature ensures a robust quantification of parameter uncertainty, thus providing reliable results.

    Network Posterior Analysis: The package offers tools to visualize and interpret network posteriors, highlighting a stable core structure and detailing individual differences and relationships.

# Installation:
The latest version of rgm can be installed from Github in R using the following command:

install.packages("devtools")
devtools::install_github("rgm/rgm")

# Usage:
For a comprehensive guide on using rgm for data analysis, refer to the package vignette and documentation:

```
library(rgm)
vignette("rgm-introduction")
```

Note: While rgm was initially designed with applications in microbiome analysis in mind, it is broadly applicable to a variety of fields that require graphical modeling of network data.
