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


## Input data into the pipeline

- Genotypes
    - plink format
    - per chromosome e.g. `<input_dir>/<prefix>1`, `<input_dir>/<prefix>2`, ..., `<input_dir>/<prefix>X`
    - Cleaned
        - MAF > ??
        - HWE > ??
        - Imputed to 1000 genomes or above
        - Info scores > ??
    - variant IDs - challenges: rsid is a positional identifier, how to handle multiple variants at one position, how to handle indels
        - chr:pos_a1_a2
            - a1 = first alphabetical
            - a2 = second alphabetical
    

## Note about variant IDs

dataset 1
A = 0.49
G = 0.51
default: 1:10000_A_G
lexical: 1:10000_A_G


dataset2
A = 0.51
G = 0.49


wrong
AA = 0
AG = 1
GG = 2

default: 1:10000_G_A
lexical: 1:10000_A_G

AA = 2
AG = 1
GG = 0



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

