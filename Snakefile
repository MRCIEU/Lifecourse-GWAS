from dotenv import dotenv_values
import glob

config = dotenv_values("config.env")

geno_input = config['genotype_input_dir']
geno_proc = config['genotype_processed_dir']
bfile_prefix = config['bfile_prefix']
resdir = config['results_dir']
phen_input = config['phenotype_input_dir']
phen_proc = config['phenotype_processed_dir']

# Get list of chr values that matches "geno_input/bfile_prefix{chr}.bed"
chr_files = glob.glob(f"{geno_input}/{bfile_prefix}*.bed")

# Extract the wildcard values from chr_files
chr = [x.split("/")[-1].split(".")[0].split(bfile_prefix)[1] for x in chr_files]

# Find the phenotypes available in the study
phen = [x.split(",")[1] for x in open("phenotype_list.csv").read().split("\n") if x]
phen.pop(0)
phen = [x for x in phen if phen_input + "/" + x + ".txt" in glob.glob(f"{phen_input}/*.txt")]

rule all:
    input:
        # expand("{geno_proc}/geno_chrs.txt", geno_proc=geno_proc),
        # expand("{geno_proc}/bfiles/vremove", geno_proc=geno_proc),
        # expand("{geno_proc}/bfiles/sremove", geno_proc=geno_proc),
        # expand("{resdir}/00/variants.txt", resdir=resdir),
        # expand("{resdir}/00/logfile", resdir=resdir),
        # expand("{resdir}/01/pcaplot.png", resdir=resdir),
        # expand("{resdir}/02/{phen}_summary.rds", resdir=resdir, phen=phen),
        # expand("{resdir}/02/organise_phenotypes.html", resdir=resdir),
        dynamic(f"{resdir}/03/{{gwas}}.fastGWA.rds")

rule check_data:
    input: 
        expand("{geno_input}/{bfile_prefix}{chr}.bed", geno_input=geno_input, bfile_prefix=bfile_prefix, chr=chr),
        expand("{geno_input}/{bfile_prefix}{chr}.bim", geno_input=geno_input, bfile_prefix=bfile_prefix, chr=chr),
        expand("{geno_input}/{bfile_prefix}{chr}.fam", geno_input=geno_input, bfile_prefix=bfile_prefix, chr=chr),
        "hello"
    output:
        expand("{geno_proc}/symlinks/{bfile_prefix}{chr}.bim", geno_proc=geno_proc, bfile_prefix=bfile_prefix, chr=chr),
        expand("{geno_proc}/geno_chrs.txt", geno_proc=geno_proc),
        expand("{geno_proc}/bfiles/vremove", geno_proc=geno_proc),
        expand("{geno_proc}/bfiles/sremove", geno_proc=geno_proc),
        expand("{resdir}/00/variants.txt", resdir=resdir),
        expand("{resdir}/00/logfile", resdir=resdir)
    shell:
        "./00-check-data.sh"

rule ancestry:
    input: 
        expand("{geno_proc}/symlinks/{bfile_prefix}{chr}.bim", geno_proc=geno_proc, bfile_prefix=bfile_prefix, chr=chr),
        expand("{geno_proc}/geno_chrs.txt", geno_proc=geno_proc),
        expand("{geno_proc}/bfiles/vremove", geno_proc=geno_proc),
        expand("{geno_proc}/bfiles/sremove", geno_proc=geno_proc),
        expand("{resdir}/00/variants.txt", resdir=resdir),
        expand("{resdir}/00/logfile", resdir=resdir)
    output:
        expand("{resdir}/01/pcaplot.png", resdir=resdir)
    shell:
        "./01-ancestry.sh"


rule phenotypes:
    input: 
        expand("{geno_proc}/{bfile_prefix}pc.txt", geno_proc=geno_proc, bfile_prefix=bfile_prefix),
        # f"{phen_input}/{{phen}}.txt",
        expand("{phen_input}/pheno_covariates.txt", phen_input=phen_input),
    output:
        # f"{resdir}/02/{{phen}}_summary.rds",
        dynamic(f"{phen_proc}/{{gwas}}.txt")
        # expand("{resdir}/02/organise_phenotypes.html", resdir=resdir),
        # expand("{resdir}/02/scores.rds", resdir=resdir),
        # expand("{resdir}/02/summary.rds", resdir=resdir)
    shell:
        "./02-phenotype-organisation.sh"


rule gwas:
    input: 
        f"{phen_proc}/{{gwas}}.txt"
    output:
        f"{resdir}/03/{{gwas}}.fastGWA.rds"
    shell:
        "./03-gwas.sh {input}"


# rule cluster:
#     input: "README.md"
#     output: dynamic("{clusterid}.something.txt")
#     run:
#         # create 5 files called 1.txt, 2.txt, 3.txt, 4.txt, 5.txt
#         for i in range(1, 6):
#             with open(f"{i}.something.txt", "w") as f:
#                 f.write("Hello, world!\n")

# rule cp:
#     input: "{clusterid}.something.txt"
#     output: "{clusterid}.something.txt2"
#     shell: "cp {input} {output}"