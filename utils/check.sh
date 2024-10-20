#!/bin/bash

source config.env


echo "Checking R packages..."
Rscript -e "renv::status()"

tf=$(mktemp)

echo "Checking plink..."
if ./bin/plink2 --version > $tf 2>&1; then
    echo "All good!"
else
    cat $tf
fi

echo "Checking flashpca..."
if ./bin/flashpca --version > $tf 2>&1; then
    echo "All good!"
else
    cat $tf
fi

echo "Checking gcta..."

./bin/gcta-1.94.1 > $tf 2>&1
if cat $tf | grep -q "Genome"; then
    echo "All good!"
else
    cat $tf
fi

echo "Checking king..."

./bin/king > $tf 2>&1
if cat $tf | grep -q "KING"; then
    echo "All good!"
else
    cat $tf
fi


