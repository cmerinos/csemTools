# CSEM: A Package for Computing Conditional Standard Errors of Measurement

The `CSEM` package provides functions to compute various methods of Conditional Standard Error of Measurement (CSEM), including:

-   **MF-CSEM**: Mollenkopf-Feldt method
-   **Thorndike-CSEM**: Based on half-test differences
-   **Binomial-CSEM**: Based on the binomial model
-   **Other Reliability-Based Methods**

Additionally, the package includes tools for: - Splitting test items into halves (`check.split`) - Computing reliability estimates (`check.alpha`, `check.angoff`) - Testing the distribution of split scores (`check.distribution`) - Visualizing results (`check.plot`)

## Installation

To install and load the package, run:

\`\`\`r \# Install from GitHub (if applicable) devtools::install_github(“cmerinos/CSEM”)

# Load the package

library(CSEM)

# Getting Started

## Example data

set.seed(123) dhalf1 \<- matrix(sample(1:5, 30, replace = TRUE), nrow = 10) dhalf2 \<- matrix(sample(1:5, 30, replace = TRUE), nrow = 10)

## Example: Computing MF-CSEM

csemMF(half1 = dhalf1, half2 = dhalf2, data = dataset, reliability.coef = 0.9, n.items = 6, min.score.item = 0, max.score.item = 4, conf.level = .95)

## Example: Computing Thorndike-CSEM

csemthorndike(half1 = dhalf1, half2 = dhalf2, n.groups = 4)

## Example: Checking Score Distributions

check.distribution(half1 = dhalf1, half2 = dhalf2, B = 2000, conf = .95)

# Documentation

?CSEM help(package = “CSEM”)

# Development

-   Author: Cesar Merino-Soto

-   Maintainer: Cesar Merino-Soto [sikayax\@yahoo.com.ar](mailto:sikayax@yahoo.com.ar){.email}

-   License: GPL-3

# Contributions

Of course, contributions always are welcome! If you find an issue, please report it or submit a pull request.
