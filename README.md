# MRcML_Family

Mendelian Randomization with cML Family Methods for Multi-Ancestry Studies

## Package

The main R package is located in the `MRcMLFamily/` subdirectory.

### Installation

```r
install.packages("devtools")
devtools::install_github("XMMMMddd/MRcML_Family/MRcMLFamily")
```

See [MRcMLFamily/README.md](MRcMLFamily/README.md) for detailed installation and usage instructions.

## Files

- `MRcMLFamily/` - R package source
  - `R/` - R source files
  - `src/` - C++ source files
  - `README.md` - Package documentation
- `MR_family.r` - Original basic MR functions
- `MRcML_family_MA.r` - Original MRcML-Family-MA implementation
- `MRcML_family_MA.cpp` - C++ acceleration module
