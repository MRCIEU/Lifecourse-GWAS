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


## Developers

To setup the environment do this:

To run tests do this:

If you wanna make changes:

1. Create an issue in github repo
2. Create a branch from that issue
3. Clone it locally
4. Make changes in that local branch
5. Push the changes to github
6. Make a pull request on github to merge to main branch

