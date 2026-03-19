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

---

## Data Format

### Beta Estimates

- **beta_hat_exp**: J x 3 matrix of exposure beta estimates (one row per SNP, three columns for three ancestries)
- **beta_hat_out**: J x 3 matrix of outcome beta estimates

### Covariance Matrices

- **beta_sigma_exp**: List of J covariance matrices (each 3x3), one for each SNP's exposure betas
- **beta_sigma_out**: List of J covariance matrices (each 3x3), one for each SNP's outcome betas

### Data Generation Helper

```r
# Helper function to generate example data
generate_example_data <- function(J = 100, seed = 42) {
  set.seed(seed)
  beta_exp <- matrix(rnorm(J * 3, sd = 0.5), nrow = J, ncol = 3)
  beta_out <- matrix(rnorm(J * 3, sd = 0.5), nrow = J, ncol = 3)
  sigma_exp <- lapply(1:J, function(j) diag(3) * runif(1, 0.01, 0.1))
  sigma_out <- lapply(1:J, function(j) diag(3) * runif(1, 0.01, 0.1))
  list(
    beta_hat_exp = beta_exp,
    beta_hat_out = beta_out,
    beta_sigma_exp = sigma_exp,
    beta_sigma_out = sigma_out
  )
}
```

---

## Function Reference

### 1. mr_family()

Basic MR family analysis with iterative estimation and variance calculation.

**Usage:**
```r
result <- mr_family(
  beta_hat_exp,
  beta_hat_out,
  beta_sigma_exp,
  beta_sigma_out,
  alpha_init = 0.0,
  max_iter = 100,
  tol = 1e-6
)
```

**Arguments:**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `beta_hat_exp` | matrix (J x 3) | - | Exposure beta estimates |
| `beta_hat_out` | matrix (J x 3) | - | Outcome beta estimates |
| `beta_sigma_exp` | list (3x3 matrices) | - | Exposure covariance matrices |
| `beta_sigma_out` | list (3x3 matrices) | - | Outcome covariance matrices |
| `alpha_init` | numeric | 0.0 | Initial value for alpha |
| `max_iter` | integer | 100 | Maximum iterations |
| `tol` | numeric | 1e-6 | Convergence tolerance |

**Returns:** Data frame with columns:
- `alpha_point`: Causal effect estimate
- `alpha_se`: Standard error
- `p_value`: Two-sided p-value

**Example:**
```r
library(MRcMLFamily)

# Generate example data
data <- generate_example_data(J = 50, seed = 123)

# Run MR analysis
result <- mr_family(
  beta_hat_exp = data$beta_hat_exp,
  beta_hat_out = data$beta_hat_out,
  beta_sigma_exp = data$beta_sigma_exp,
  beta_sigma_out = data$beta_sigma_out
)

# View results
print(result)
#   alpha_point   alpha_se   p_value
# 1    0.2841    0.3425    0.4068
```

---

### 2. MRcML_family()

Main analysis function with BIC-based Bayesian Model Averaging across (a, b) parameter grid.

**Usage:**
```r
result <- MRcML_family(
  a_values,
  b_values,
  beta_hat_exp,
  beta_hat_out,
  beta_sigma_exp,
  beta_sigma_out,
  n = 1000,
  p = 1,
  q = 1/3
)
```

**Arguments:**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `a_values` | integer vector | 1:n_snps | Grid of a parameters (horizontal pleiotropy) |
| `b_values` | integer vector | 0:n_snps | Grid of b parameters (vertical pleiotropy) |
| `beta_hat_exp` | matrix (J x 3) | - | Exposure beta estimates |
| `beta_hat_out` | matrix (J x 3) | - | Outcome beta estimates |
| `beta_sigma_exp` | list (3x3 matrices) | - | Exposure covariance matrices |
| `beta_sigma_out` | list (3x3 matrices) | - | Outcome covariance matrices |
| `n` | integer | 1000 | Sample size for BIC penalty |
| `p` | numeric | 1 | Number of exposures |
| `q` | numeric | 1/3 | Ratio parameter |

**Returns:** List with components:
- `alpha_bma`: BMA causal effect estimate
- `se_bma`: BMA standard error
- `results_with_weights`: Full grid results data frame

**Example:**
```r
library(MRcMLFamily)

# Generate example data
data <- generate_example_data(J = 100, seed = 456)
n_snps <- nrow(data$beta_hat_exp)

# Run MRcML-Family analysis
result <- MRcML_family(
  a_values = 1:n_snps,
  b_values = 0:n_snps,
  beta_hat_exp = data$beta_hat_exp,
  beta_hat_out = data$beta_hat_out,
  beta_sigma_exp = data$beta_sigma_exp,
  beta_sigma_out = data$beta_sigma_out,
  n = 1000,
  p = 1,
  q = 1/3
)

# Access results
cat("BMA Alpha:", result$alpha_bma, "\n")
cat("BMA SE:", result$se_bma, "\n")

# View grid results
head(result$results_with_weights)
```

---

### 3. cml_family_new()

Alias for `MRcML_family()` with identical functionality.

**Usage:**
```r
result <- cml_family_new(
  a_values,
  b_values,
  beta_hat_exp,
  beta_hat_out,
  beta_sigma_exp,
  beta_sigma_out,
  n = 1000,
  p = 1,
  q = 1/3
)
```

**Example:**
```r
library(MRcMLFamily)

# Same as MRcML_family but with different name
data <- generate_example_data(J = 80, seed = 789)
n_snps <- nrow(data$beta_hat_exp)

result <- cml_family_new(
  a_values = 1:n_snps,
  b_values = 0:n_snps,
  beta_hat_exp = data$beta_hat_exp,
  beta_hat_out = data$beta_hat_out,
  beta_sigma_exp = data$beta_sigma_exp,
  beta_sigma_out = data$beta_sigma_out
)
```

---

### 4. run_cml_family()

Core iterative algorithm for cML family estimation with fixed (a, b).

**Usage:**
```r
result <- run_cml_family(
  beta_hat_exp,
  beta_hat_out,
  beta_sigma_exp,
  beta_sigma_out,
  a = 10,
  b = 10,
  alpha_init = 0.0,
  max_iter = 100,
  tol = 1e-6
)
```

**Arguments:**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `beta_hat_exp` | matrix (J x 3) | - | Exposure beta estimates |
| `beta_hat_out` | matrix (J x 3) | - | Outcome beta estimates |
| `beta_sigma_exp` | list (3x3 matrices) | - | Exposure covariance matrices |
| `beta_sigma_out` | list (3x3 matrices) | - | Outcome covariance matrices |
| `a` | integer | 10 | Number of top SNPs for exposure |
| `b` | integer | 10 | Number of top SNPs for outcome |
| `alpha_init` | numeric | 0.0 | Initial alpha value |
| `max_iter` | integer | 100 | Maximum iterations |
| `tol` | numeric | 1e-6 | Convergence tolerance |

**Returns:** List with:
- `alpha`: Estimated causal effect
- `gamma`: Estimated gamma values
- `S_f`: Selected exposure SNP indices
- `S_t`: Selected outcome SNP indices

**Example:**
```r
library(MRcMLFamily)

# Generate example data
data <- generate_example_data(J = 50, seed = 111)

# Run with specific (a, b)
result <- run_cml_family(
  beta_hat_exp = data$beta_hat_exp,
  beta_hat_out = data$beta_hat_out,
  beta_sigma_exp = data$beta_sigma_exp,
  beta_sigma_out = data$beta_sigma_out,
  a = 5,
  b = 10,
  alpha_init = 0.0
)

# View results
cat("Alpha estimate:", result$alpha, "\n")
cat("Selected exposure SNPs (S_f):", result$S_f, "\n")
cat("Selected outcome SNPs (S_t):", result$S_t, "\n")
```

---

## Understanding Parameters

- **a (horizontal pleiotropy)**: Number of top SNPs to select based on exposure statistics. Controls the degree of horizontal pleiotropy modeling.
- **b (vertical pleiotropy)**: Number of top SNPs to select based on outcome statistics. Controls vertical pleiotropy modeling.
- **n**: Sample size for BIC penalty calculation
- **p**: Number of exposures
- **q**: Ratio parameter for model complexity penalty (default 1/3)

### Parameter Selection Guide

- For small SNP sets (J < 50): Use `a_values = 1:J`, `b_values = 0:J`
- For large SNP sets (J > 100): Use smaller grid, e.g., `a_values = c(5, 10, 20, 50)`, `b_values = c(0, 5, 10, 20)`
- Larger a/b values increase model complexity; use BIC to select optimal

---

## Output Interpretation

### Basic MR (`mr_family()`)

Returns a data frame:
- `alpha_point`: Point estimate of causal effect
- `alpha_se`: Standard error
- `p_value`: Two-sided p-value (from normal distribution)

### Advanced MRcML-Family (`MRcML_family()`)

Returns a list:
- `alpha_bma`: Bayesian Model Averaging estimate
- `se_bma`: BMA standard error (includes model uncertainty)
- `results_with_weights`: Full data frame with columns:
  - `a`, `b`: Parameter values
  - `alpha_hat`: Point estimate for that (a,b)
  - `variance`, `std_error`: Variance estimates
  - `p_value`: Two-sided p-value
  - `bic`: Bayesian Information Criterion
  - `weight`: BIC weight for BMA

---

## Complete Workflow Example

```r
library(MRcMLFamily)

# 1. Prepare your data (replace with real GWAS summary statistics)
data <- generate_example_data(J = 100, seed = 42)
n_snps <- nrow(data$beta_hat_exp)

# 2. Quick analysis with mr_family()
quick_result <- mr_family(
  beta_hat_exp = data$beta_hat_exp,
  beta_hat_out = data$beta_hat_out,
  beta_sigma_exp = data$beta_sigma_exp,
  beta_sigma_out = data$beta_sigma_out
)
print(quick_result)

# 3. Comprehensive analysis with MRcML_family()
full_result <- MRcML_family(
  a_values = 1:n_snps,
  b_values = 0:n_snps,
  beta_hat_exp = data$beta_hat_exp,
  beta_hat_out = data$beta_hat_out,
  beta_sigma_exp = data$beta_sigma_exp,
  beta_sigma_out = data$beta_sigma_out,
  n = 1000,
  p = 1,
  q = 1/3
)

# 4. Extract results
cat("\n=== MRcML-Family Results ===\n")
cat("BMA Alpha:", full_result$alpha_bma, "\n")
cat("BMA SE:", full_result$se_bma, "\n")

# 5. Find best model by BIC weight
best_model <- full_result$results_with_weights[
  which.max(full_result$results_with_weights$weight),
]
cat("\nBest model (a, b):", best_model$a, best_model$b, "\n")
cat("BIC weight:", best_model$weight, "\n")
```

---

## Methods

The package implements:

1. **Iterative Optimization**: EM-style algorithm for joint estimation of alpha and gamma
2. **Model Selection**: BIC-based selection across (a, b) parameter grid
3. **Robust Estimation**: Multi-start strategies with BIC-weighted aggregation
4. **Variance Estimation**: Fisher information-based standard errors
5. **Bayesian Model Averaging**: Weighted combination across candidate models

## License

MIT
