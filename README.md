# Lifecourse-GWAS

- config file in json
- if we use simple bash scripting then include jq binary in the /bin folder, but will need to detect OS and download the correct version
- alternative to use snakemake / R/targets / only R files
- keep data and outputs separate
- need example / test data


## Setup environment

In R:
```r
renv::restore()
```

## Run

1. Phenotype organisation

```r
Rscript resources/render.r resources/phenotypes/organise_phenotypes.rmd
```
