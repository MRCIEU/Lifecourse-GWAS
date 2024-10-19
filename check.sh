#!/bin/bash

source config.env

Rscript -e "renv::status()"

./bin/plink2 --version
./bin/flashpca --version
./bin/gcta-1.94.1
./bin/king
