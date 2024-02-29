#!/bin/bash

# Get model
grep -riL "use linear" *log
grep "use linear" *log | cut -d ":" -f 1

# Get sample size
grep -m 1 "overlapping" *log | cut -d " " -f 1

# Get h2

# Get runtime

# Number of SNPs

