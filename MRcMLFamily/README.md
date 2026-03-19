# MRcMLFamily

Mendelian Randomization with cML Family Methods for Multi-Ancestry Studies

## Overview

`MRcMLFamily` is an R package implementing Mendelian Randomization (MR) methods using collaborative IV (cML) family approaches. It is designed for handling horizontal pleiotropy in multi-ancestry genetic studies.

## Installation

### From GitHub

```r
# Install devtools if needed
install.packages("devtools")

# Install MRcMLFamily from GitHub
devtools::install_github("XMMMMddd/MRcML_Family/MRcMLFamily")
```

### From Local Source

```r
# Clone the repository first
# git clone https://github.com/XMMMMddd/MRcML_Family.git

# Install from local source
devtools::install("/path/to/MRcML_Family/MRcMLFamily")
```

### Build from Source (Terminal)

```bash
# Navigate to package directory
cd MRcML_Family/MRcMLFamily

# Build package
R CMD build .

# Install
R CMD INSTALL MRcMLFamily_1.0.0.tar.gz
```

## Requirements

- R (>= 4.0)
- Rcpp (>= 1.0.10)
- RcppArmadillo

These dependencies will be automatically installed when you install the package.

## Quick Start

```r
library(MRcMLFamily)

# Prepare your data (see below for data format)
# beta_exp, beta_out: J x 3 matrices of beta estimates
# sigma_exp, sigma_out: lists of 3x3 covariance matrices

# Basic MR analysis
result <- mr_family(
  beta_hat_exp = beta_exp,
  beta_hat_out = beta_out,
  beta_sigma_exp = sigma_exp,
  beta_sigma_out = sigma_out
)

# View results
print(result)
```

## Data Format

### Beta Estimates

- **beta_hat_exp**: J x 3 matrix of exposure beta estimates (one row per SNP, three columns for three ancestries)
- **beta_hat_out**: J x 3 matrix of outcome beta estimates

### Covariance Matrices

- **beta_sigma_exp**: List of J covariance matrices (each 3x3), one for each SNP's exposure betas
- **beta_sigma_out**: List of J covariance matrices (each 3x3), one for each SNP's outcome betas

### Example Data Generation

```r
# Number of SNPs
J <- 100
n_snps <- J

# Generate beta estimates
set.seed(42)
beta_exp <- matrix(rnorm(J * 3, sd = 0.5), nrow = J, ncol = 3)
beta_out <- matrix(rnorm(J * 3, sd = 0.5), nrow = J, ncol = 3)

# Generate covariance matrices
sigma_exp <- lapply(1:J, function(j) {
  diag(3) * runif(1, 0.01, 0.1)
})
sigma_out <- lapply(1:J, function(j) {
  diag(3) * runif(1, 0.01, 0.1)
})
```

## Main Functions

### Basic MR Analysis

| Function | Description |
|----------|-------------|
| `mr_family()` | Basic MR family analysis with variance estimation |

### Advanced MRcML-Family Analysis

| Function | Description |
|----------|-------------|
| `MRcML_family()` | Main analysis with BIC-based model averaging |
| `cml_family_new()` | Alias for `MRcML_family()` |
| `run_cml_family()` | Core iterative algorithm |

### Utility Functions

| Function | Description |
|----------|-------------|
| `run_cml_family_internal()` | Internal iterative algorithm |
| `calculate_alpha_variance_internal_two()` | Variance estimation |
| `calculate_bic()` | BIC calculation |
| `tune_cml_family_grid_fast()` | Grid search over (a,b) parameters |
| `calculate_bic_weighted_estimates_fast()` | BIC-based Bayesian Model Averaging |

## Usage Examples

### Example 1: Basic MR Analysis

```r
library(MRcMLFamily)

# Generate example data
set.seed(123)
J <- 50
beta_exp <- matrix(rnorm(J * 3, sd = 0.3), nrow = J, ncol = 3)
beta_out <- matrix(rnorm(J * 3, sd = 0.3), nrow = J, ncol = 3)
sigma_exp <- lapply(1:J, function(j) diag(3) * 0.05)
sigma_out <- lapply(1:J, function(j) diag(3) * 0.05)

# Run basic MR analysis
result <- mr_family(
  beta_hat_exp = beta_exp,
  beta_hat_out = beta_out,
  beta_sigma_exp = sigma_exp,
  beta_sigma_out = sigma_out
)

# Results
cat("Alpha estimate:", result$alpha_point, "\n")
cat("Standard error:", result$alpha_se, "\n")
cat("P-value:", result$p_value, "\n")
```

### Example 2: Full MRcML-Family Analysis

```r
library(MRcMLFamily)

# Generate example data
set.seed(456)
J <- 100
beta_exp <- matrix(rnorm(J * 3, sd = 0.4), nrow = J, ncol = 3)
beta_out <- matrix(rnorm(J * 3, sd = 0.4), nrow = J, ncol = 3)
sigma_exp <- lapply(1:J, function(j) diag(3) * 0.02)
sigma_out <- lapply(1:J, function(j) diag(3) * 0.02)

# Run MRcML-Family with BIC model averaging
n_snps <- nrow(beta_exp)
result <- MRcML_family(
  a_values = 1:n_snps,
  b_values = 0:n_snps,
  beta_hat_exp = beta_exp,
  beta_hat_out = beta_out,
  beta_sigma_exp = sigma_exp,
  beta_sigma_out = sigma_out,
  n = 1000,
  p = 1,
  q = 1/3
)

# Access results
cat("BMA Alpha estimate:", result$alpha_bma, "\n")
cat("BMA Standard error:", result$se_bma, "\n")

# View full grid results
head(result$results_with_weights)
```

### Example 3: Custom (a, b) Parameter Grid

```r
library(MRcMLFamily)

# Use a smaller grid for faster computation
result <- MRcML_family(
  a_values = c(5, 10, 15, 20),
  b_values = c(0, 5, 10, 15),
  beta_hat_exp = beta_exp,
  beta_hat_out = beta_out,
  beta_sigma_exp = sigma_exp,
  beta_sigma_out = sigma_out,
  n = 1000,
  p = 1,
  q = 1/3
)
```

## Understanding Parameters

- **a (horizontal pleiotropy)**: Number of top SNPs to select based on exposure statistics
- **b (vertical pleiotropy)**: Number of top SNPs to select based to outcome statistics
- **n**: Sample size for BIC penalty calculation
- **p**: Number of exposures
- **q**: Ratio parameter for model complexity penalty

## Output Interpretation

The main analysis returns a list with:

- **alpha_bma**: Bayesian Model Averaging (BMA) estimate of the causal effect
- **se_bma**: BMA standard error accounting for model uncertainty
- **results_with_weights**: Dataframe with results from all (a,b) combinations, including:
  - `a`, `b`: Parameter values
  - `alpha_hat`: Point estimate for that (a,b) combination
  - `variance`, `std_error`: Variance and standard error
  - `p_value`: Two-sided p-value
  - `bic`: Bayesian Information Criterion
  - `weight`: BIC weight for BMA aggregation

## Methods

The package implements:

1. **Iterative Optimization**: EM-style algorithm for joint estimation of alpha and gamma
2. **Model Selection**: BIC-based selection across (a, b) parameter grid
3. **Robust Estimation**: Multi-start strategies with BIC-weighted aggregation
4. **Variance Estimation**: Fisher information-based standard errors
5. **Bayesian Model Averaging**: Weighted combination across candidate models

## License

MIT
