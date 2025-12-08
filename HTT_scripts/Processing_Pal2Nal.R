# Set working directory
setwd("~/Downloads/")

library(seqinr)

##splitting RVT to individual fastas
# Load the full FASTA file
seqs <- read.fasta("rvt_full_nt.fasta", seqtype = "DNA", as.string = TRUE)

# Create output folder
dir.create("split_fastas", showWarnings = FALSE)

# Write each sequence to its own FASTA file
for (name in names(seqs)) {
  file_name <- paste0("split_fastas/", name, ".fasta")
  write.fasta(sequences = list(seqs[[name]]),
              names = name,
              file.out = file_name)
}

##create TE-TE seq files from tsv
# Load TE pair table (ignore 3rd column)
te_pairs <- read.table("species-unique_filtered_pim_pairs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(te_pairs) <- c("TE1", "TE2", "Ignore")

# Create output directory
dir.create("pair_fastas", showWarnings = FALSE)

# Loop through each pair
for (i in 1:nrow(te_pairs)) {
  te1 <- te_pairs$TE1[i]
  te2 <- te_pairs$TE2[i]
  
  # Construct file paths
  file1 <- paste0("split_fastas/", te1, ".fasta")
  file2 <- paste0("split_fastas/", te2, ".fasta")
  
  # Read sequences
  seq1 <- read.fasta(file = file1, seqtype = "DNA", as.string = TRUE)
  seq2 <- read.fasta(file = file2, seqtype = "DNA", as.string = TRUE)
  
  # Output filename based on pair name
  out_file <- paste0("pair_fastas/", te1, "_vs_", te2, ".fasta")
  
  # Write combined file
  write.fasta(sequences = c(seq1, seq2),
              names = c(names(seq1), names(seq2)),
              file.out = out_file)
}

#KaKs

library(seqinr)

# Set working directory
setwd("~/Downloads/codon_alignments_fasta/")

# Load the existing TSV file
pair_data <- read.table("species-unique_filtered_pim_pairs.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Add new columns for Ka, Ks, and Ka/Ks if not already present
if (!all(c("Ka", "Ks", "Ka_Ks_ratio") %in% colnames(pair_data))) {
  pair_data$Ka <- NA
  pair_data$Ks <- NA
  pair_data$Ka_Ks_ratio <- NA
}

# List all .phy codon alignment files
phy_files <- list.files(pattern = "\\.phy$")

# Function to extract pair name from filename (e.g., CalAnn_X_vs_Y_codon_alignment.phy)
extract_pair <- function(filename) {
  gsub("_codon_alignment\\.phy$", "", filename)
}

# Loop through each file
for (file in phy_files) {
  cat("Processing:", file, "\n")
  try({
    aln <- read.alignment(file, format = "fasta")
    result <- kaks(aln)
    
    ka <- result$ka
    ks <- result$ks
    ratio <- result$ka.ks
    
    pair_name <- extract_pair(file)
    
    # Match and fill into the pair_data table (match Row vs Column combinations)
    pair_data_match <- grepl(pair_name, paste(pair_data$Row, pair_data$Column, sep = "_vs_")) |
      grepl(pair_name, paste(pair_data$Column, pair_data$Row, sep = "_vs_"))
    
    if (any(pair_data_match)) {
      pair_data$Ka[pair_data_match] <- ka
      pair_data$Ks[pair_data_match] <- ks
      pair_data$Ka_Ks_ratio[pair_data_match] <- ratio
    } else {
      cat("Warning: No matching pair found for", pair_name, "\n")
    }
  }, silent = TRUE)
}

# Save updated table
write.table(pair_data, file = "species-unique_filtered_pim_pairs_with_kaks.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
