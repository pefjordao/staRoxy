# staRoxy `0.1.0` <img src="staRoxy/vignettes/staRoxy_logo.svg" align="right" width="180">

**staRoxy** is a dedicated R package designed to streamline oxylipidomics abundance data analysis. It provides a reproducible framework for data cleaning, normalization, statistical modeling, and integrated visualization. While originally developed for oxylipin profiling (in press), it is adaptable for other omics datasets, including transcriptomics, proteomics, and metabolomics. 

---

## Installation

**staRoxy** is currently not available on the Comprehensive R Archive Network (CRAN). The development version can be installed via **BiocManager**, which will also install the required dependencies:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("pefjordao/staRoxy", build_vignettes = TRUE)
```

## Documentation

For a complete step-by-step guide, including data requirements and statistical assumptions, please refer to the official vignette:

<https://pefjordao.github.io/staRoxy/>

## Author/Support

Pedro Henrique F. Jordão

pedro.hf.jordao@usp.br

<https://github.com/pefjordao/staRoxy/issues>

## Funding

The development of the staRoxy package has been partially supported by:

- São Paulo Research Foundation (FAPESP): grant number 2025/00242-1.
- FAPESP/CEPID-Redoxoma: grant number 2013/07937-8.
