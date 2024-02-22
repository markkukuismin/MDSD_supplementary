# MDSD_supplementary

Supplementary materials for the article 'Network hub detection using the entire solution path information'.

See also [https://github.com/markkukuismin/MDSD](https://github.com/markkukuismin/MDSD).

## Dependencies

```r
install.packages(c("igraph", "ggplot2", "huge", "hglasso", "tidyverse", "gridExtra"))
```

Package `space` package is removed from the CRAN repository (last visit 20/2/2024) and must be installed using the tar.gz file.

Packages `DESeq2` and `genefilter` are available on Bioconductor,

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "genefilter"))
```

## Main article Figures

Use these scripts to reproduce Figures of the article,

* `manuscript_figures/Fig1.R`: Main article Figure 1.
* `manuscript_figures/Supp_Fig2_and_Fig3.R`: Main article Figure 2.
* `real_data_examples/Maize_ligule/Maize_ligule_hglasso.R`: Main article Figure 3.

## Simulations

Use these scripts to reproduce the simulations results reported in the main article and Supplementary materials,

* `simulations/simulations_lossy_screening.R`: Lossy screening rule of the sample correlation matrix simulations.
* `simulations/simulations_hglasso.R`: Hub glasso simulations.
* `simulations/simulations_local_hub_screening.R`: Firouzi & Hero (2013) hub screening method simulations.
* `simulations/simulations_space.R`: space simulations.

## Real data example

Use the script to reproduce the real data example reported in the main article. The data set is publicly available at [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61333](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61333)

* `Maize_ligule_hglasso.R`: Derive a gene co-expression network to detect hub genes of Maize ligule data.
