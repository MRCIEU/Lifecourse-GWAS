# Example input data

The data provided here is randomly generated and will not necessarily reflect realistic values. They are provided for the purposes of unit testing and to illustrate the formats of the various files expected by the pipeline.

Notes on particular files:

- `config.env`: Note that an extra `static_covariate` has been added (`batch`) for illustrative purposes. This `batch` variable is hence included in the `static_covariates.txt` file.
- `geno_input.txt`: Lists the locations of the per-chromosome `bgen` and `sample` files.
- `phen_input/static_covariates.txt`: Contains the static covariates for the individuals in the study. The columns correspond to the columns noted in the `static_covariates` field in `config.env`. Note that this is not time-varying (i.e. age is not required here).
- `phen_input/`: This directory is where the pipeline expects to find all the phenotype and covariate files. The pipeline will look for the files listed in the [phenotype_list.csv](https://github.com/MRCIEU/Lifecourse-GWAS/blob/main/phenotype_list.csv) file in this directory. The `config.env` file specifies the location of this directory in the `phenotype_input_dir` field.
- `phen_input/bmi.txt`: Contains the BMI values for the individuals in the study, along with the age at measurement. Note that individuals can have multiple measurements at different ages, though this is not required.
- `phen_input/shbg.txt`: Contains the SHBG values for the individuals in the study, along with the age at measurement. As specified in [phenotype_list.csv](https://github.com/MRCIEU/Lifecourse-GWAS/blob/main/phenotype_list.csv), this phenotype has specific covariates required, and correspondingly the `hormone_med` column is present in the `phen_input/shbg.txt` phenotype file.


