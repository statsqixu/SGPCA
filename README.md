# SGPCA: Sparse Group Principal Component Analysis

## Overview

SGPCA is an R package that implements Sparse Group Principal Component Analysis using a double thresholding algorithm. This method combines both individual and group-level sparsity in principal component analysis, making it particularly effective for analyzing high-dimensional data with natural group structures (e.g., gene expression data across cell types).

## Features

- **Double Thresholding Algorithm**: Implements both individual and group-level sparsity
- **Flexible PC Estimation**: Support for single or multiple principal components
- **Parameter Tuning**: Includes resampling-based automatic parameter selection
- **Interactive Mode**: Offers manual parameter selection with visualization
- **Built-in Simulators**: Three different simulation settings for testing and validation

## Installation

```R
install.packages("devtools")
devtools::install_github("statsqixu/SGPCA")
```

## Example

```R
library(SGPCA)

# Generate data
data <- simulator1(n = 100, G = 300, C = C, seed = seed)
X <- data$X
group_label <- data$group_label
signal_indices <- data$signal_indices
pc1 <- data$pc1

# Estimate the first PC via SGPCA and select tuning parameters via resampling
SGPCA_results <- SGPCA.rs(X, group_label, J = 1, B = 20, rho = 0.5,
                                    tau_range = 10 ^ seq(-2, 3, length.out = 20),
                                    eta_range = 10 ^ seq(-2, 3, length.out = 20),
                                    max_iter = 20, tol = 1e-5, verbose = FALSE, mode = "auto")

```

## References

[To be updated]