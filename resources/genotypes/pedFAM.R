#!/usr/bin/env RScript
#' 1. This script will generate the FAM for fastFAM using provided pedigree information
#'   (like the ukb_inferredRelationship.txt by UKB cohort). 
#' 2. Two input files are needed: 
#'    a. A .fam file from the PLINK format genotype. It contains the individuals' FIDs and IIDs and some other information. 
#'        see (https://www.cog-genomics.org/plink/1.9/formats#fam) for details.
#'    b. A inferredRelationship file. We will use the UKB provided one. 
#'        It should contains: "IID1", "IID2", and a "relationship" columns. The "relationship" column should contain
#'        the following elements (these are used by UK Biobank to indicate the inferred relationship between 
#'        two individuals, see https://www.biorxiv.org/content/early/2017/07/20/166298) :
#'          "MZtwin", 
#'          "parentOffspring",
#'          "fullSib", 
#'          "second",
#'          "third"
#'          
#'  Author: Longda Jiang, longda.jiang@uq.edu.au
#'  Date: June 28, 2018
#'  Version: 1.0.1
#'  Modified by Gibran Hemani to accept King output (2024-02-07)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3){
    message("Input: path of fam file")
    message("Input: path of pedigree information")
    message("Output: path of output without extension")
    stop("The FAM from pedigree needs 3 parameters")
}

path_to_fam <- as.character(args[1])
path_to_rel <- as.character(args[2])
output_file <- as.character(args[3])
message("Transform the pedigree information...")
message("From: ",  path_to_rel)
message("To: ", output_file)
message("Align order: ", path_to_fam)

if(!file.exists(path_to_fam)){
    stop("FAM file didn't exist")
}

if(!file.exists(path_to_rel)){
    stop("Relationship file didn't exist")
}


###############
# Load the input data
###############
# the .fam file
fam_file <- read.table(file=path_to_fam, stringsAsFactors = F, head = F )
fam_file <- fam_file[, 1:2]
names(fam_file) <- c("FID", "IID")
fam_file$order <- 1:nrow(fam_file)

rel <- read.table(file=path_to_rel, head=F, stringsAsFactors = F) 
# the rel file we only need the 1, 2, and last columns:
rel <- rel[, c(1, 2, ncol(rel))]
names(rel) <- c("IID_1", "IID_2", "relationship")

# check relation labels
relation_labels = table(rel$relationship)
print(relation_labels)
setting_labels = c("Dup/MZ", "PO", "FS", "2nd", "3rd")
rel <- subset(rel, relationship %in% setting_labels)

rel2 <- rel[,c(2,1,3)]
names(rel2) <- c("IID_1", "IID_2", "relationship")
rel <- rbind(rel, rel2)

###############
# format the data
###############
rel <- merge(rel, fam_file[,c("IID", "order")], by.x="IID_1", by.y="IID")
rel <- merge(rel, fam_file[,c("IID", "order")], by.x="IID_2", by.y="IID")
rel <- rel[, c("IID_1", "IID_2", "relationship", "order.x", "order.y")]

# The order of IID_1 is important: 
# The order of IID_1 should always be larger than the order of IID_2. 
#  (as the sparse FAM is saved as a lower triangular matrix)
rel <- rel[order(rel$order.x, rel$order.y), ]
sum(rel$order.x >= rel$order.y)
rel <- rel[rel$order.x > rel$order.y, ]




###############
# generate the FAM
###############

rel$coef <- 0 

rel[rel$relationship == "Dup/MZ", "coef"] <- 1
rel[rel$relationship == "PO", "coef"] <- 0.5
rel[rel$relationship == "FS", "coef"] <- 0.5
rel[rel$relationship == "2nd", "coef"] <- 0.25
rel[rel$relationship == "3rd", "coef"] <- 0.125

# the diagonal elements of FAM is simply 1
rel_part2 <- fam_file[,c("IID", "IID")]
rel_part2$order.x <- 1:nrow(rel_part2)
rel_part2$order.y <- 1:nrow(rel_part2)
names(rel_part2) <- c("IID_1", "IID_2", "order.x", "order.y")
rel_part2$coef <- 1

####
FAM_sp <- rbind(rel[, c("IID_1", "IID_2", "order.x", "order.y", "coef")], rel_part2)
FAM_sp <- FAM_sp[order(FAM_sp$order.x, FAM_sp$order.y), ]

### the actual order starts from 0. need to minus 1. 
FAM_sp$order.x <- FAM_sp$order.x - 1
FAM_sp$order.y <- FAM_sp$order.y - 1

write.table(fam_file, file=paste0(output_file, ".grm.id"), quote=F, row.names = F, col.names = F)
write.table(FAM_sp[,c("order.x", "order.y", "coef")], file=paste0(output_file, ".grm.sp"), quote=F, row.names = F, col.names = F)

message("Sparse FAM generated: ", output_file, ".")

####  Done

# Greedy remove related individuals

message("Generating list of unrelated individuals...")

t1 <- Sys.time()
r <- subset(rel, coef > 0.05)
relids <- unique(c(r$IID_1, r$IID_2))

i <- 1
n <- nrow(r)
while(nrow(r) > 0) {
    # message(nrow(r)/n)
    if((nrow(r) / n) <= ((10-i) / 10)) {
        message("Progress: ", round(i / 10 * 100, 2), "%")
        i <- i + 1
    }
    a <- table(c(r$IID_1, r$IID_2))
    id <- names(a)[which.max(a)][1]
    relids <- relids[relids != id]
    r <- subset(r, !(IID_1 %in% id | IID_2 %in% id))
}
Sys.time() - t1

allrels <- unique(c(rel$IID_1, rel$IID_2))
to_remove <- allrels[!allrels %in% relids]

message("Individuals to remove: ", length(to_remove))
fam_file_unrelated <- subset(fam_file, !(IID %in% to_remove))
message("Unrelated individuals: ", nrow(fam_file_unrelated))

write.table(fam_file_unrelated[,1:2], file=paste0(output_file, ".unrelated"), quote=F, row.names = F, col.names = F)
