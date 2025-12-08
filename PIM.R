####PIM####
##important first step to generate pairwise alignments of all TEs to each other AND then filter to keep only those with 75% or higher identitiy

#run this first
#clustalo -i X.clw --percent-id --distmat-out=pim.txt --full --force

# Load necessary libraries
install.packages("reshape2")
library(reshape2)

# Set working directory
setwd("~/Documents/CR1/PIM/")

### DO FIRST ###
#remove first cell value generated in the PIM
#transpose the headers so both col and row have the same headers (around the PIM itself)
#save as .tsv file

# Read in the PIM matrix
pim_matrix <- read.table("pim.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Convert to matrix to be safe
pim_matrix <- as.matrix(pim_matrix)

# Melt into long format (this is key!)
long_pim <- reshape2::melt(pim_matrix, varnames = c("Row", "Column"), value.name = "PIM")

# Filter for PIM > 75 & < 100, and remove self-hits
filtered_pairs <- subset(long_pim, PIM > 75 & PIM < 100 & Row != Column)

# Create unique PairID (e.g., A_B or B_A becomes the same)
filtered_pairs$PairID <- apply(filtered_pairs[, c("Row", "Column")], 1, function(x) {
  paste(sort(x), collapse = "_")
})

# Remove duplicated pairs (A_B vs B_A)
unique_pairs <- filtered_pairs[!duplicated(filtered_pairs$PairID), ]

# Drop the helper column
unique_pairs$PairID <- NULL

# Write to file
write.table(unique_pairs, "unique_filtered_pim_pairs.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Extract species names from sequence IDs (before the first "_")
get_species <- function(x) sub("_.*", "", x)

# Add species columns
filtered_pairs$RowSpecies <- get_species(filtered_pairs$Row)
filtered_pairs$ColSpecies <- get_species(filtered_pairs$Column)

# Filter out same-species pairs
filtered_pairs <- filtered_pairs[filtered_pairs$RowSpecies != filtered_pairs$ColSpecies, ]

# Continue deduplication
filtered_pairs$PairID <- apply(filtered_pairs[, c("Row", "Column")], 1, function(x) {
  paste(sort(x), collapse = "_")
})
unique_pairs <- filtered_pairs[!duplicated(filtered_pairs$PairID), ]

# Drop helper columns
unique_pairs$PairID <- NULL
unique_pairs$RowSpecies <- NULL
unique_pairs$ColSpecies <- NULL

# Write to file
write.table(unique_pairs, "species-unique_filtered_pim_pairs.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

### SEQKIT ###
#get header IDs from species-unique_filtered_pim_pairs.tsv for both columns
#seqkit grep -f rvt-headers.txt rvt_nt.fasta > rvt_nt_75.fasta
#seqkit grep -f rvt-headers.txt rvt_aa.fasta > rvt_aa_75.fasta
