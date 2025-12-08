####KAKS FILTER####
##The all.dNdS,txt depends on the output of 05-coreGenedS.R from HeloiseMuller/HTvertebrates

library(dplyr)
setwd(path.expand("~/Documents/CR1/dS/dNdS/"))

# read table from text file
df <- read.table("all.dNdS.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# extract just the species names
df <- df %>%
  mutate(
    query_species   = sub(".*[\\|][+-]:(.*)$", "\\1", query),
    subject_species = sub(".*[\\|][+-]:(.*)$", "\\1", subject),
    # symmetric species pair (ignore scaffold IDs)
    species_pair    = apply(cbind(query_species, subject_species), 1, function(x) paste(sort(x), collapse = "_vs_"))
  )

# filtering
df_filtered <- df %>%
  filter(dS < 9,
         alnLength >= 100,
         dN != 0,
         !is.na(dS))

# quantiles per *species pair* (all scaffolds collapsed)
quantiles <- df_filtered %>%
  group_by(species_pair) %>%
  summarise(
    count       = n(),
    dS_q0.5perc = quantile(dS, probs = 0.005),
    dS_median   = quantile(dS, probs = 0.5),
    .groups = "drop"
  )

# save summary table
write.table(quantiles, "species_pair_quantiles.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# save filtered dataframe to a text file
write.table(df_filtered,
            file = "filtered_alignments.txt",
            sep = "\t",
            row.names = FALSE,at a notably higher rate
            quote = FALSE)
