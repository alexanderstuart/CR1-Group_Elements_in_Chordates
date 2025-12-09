
library(ape)
library(tidyverse)
library(ggplot)

group <- "CR1_BEAST_aligned"

# name_index_dir <- "../_general_data/Transposons/Colab/Alex/CR1_project/Working/RVT/final_alignment_parsed_index"
name_index_dir <- "../_general_data/Transposons/Colab/Alex/CR1_project/Working/RVT/BEAST_alignment_parsed_index"

# clade_index_dir <- "data/Phylogeny/CR1_project/clade_labels"
clade_index_dir <- "data/Phylogeny/CR1_project/BEAST_clade_labels"


data_dir <- paste0("data/parse/Cluster_join/", group, "/top_edges")

# --- Load data ---------------------------------------------------------------
df <- read_delim(data_dir, col_names = F) %>% 
      setNames(c("from", "to", "bit", "id"))

# df <- read_delim("data/parse/Cluster_join/CR1_final_aligned/external_top_edges", col_names = F) %>% 
#       setNames(c("from", "to", "bit", "id"))

# df <- read_delim("data/parse/Cluster_join/CR1_final_aligned/all_edges", col_names = F) %>% 
#       setNames(c("from", "to", "bit", "id"))

name_index <- read_delim(name_index_dir, col_names = F) %>% 
              setNames(c("new", "old"))

clade_index <- read_delim(clade_index_dir, col_names = F) %>% 
              setNames(c("node", "clade", "clade_number"))

tr <- read.tree("data/Phylogeny/CR1_project/final_rvt.nwk")

# --- merge data -------------------------------------------
df_parsed <- df %>% mutate(from_superfamily = str_split_i(from, "_", 2)) %>% 
                    mutate(to_superfamily = str_split_i(to, "_", 2)) %>%
                    mutate(from_superfamily = str_replace(from_superfamily, "-Babar", "1")) %>%
                    mutate(to_superfamily = str_replace(to_superfamily, "-Babar", "1"))

df_parsed <- df_parsed %>% left_join(name_index, by = c("from" = "old" )) %>%
                           left_join(name_index, by = c("to" = "old" )) %>%        
                           select(-from, -to) %>%
                           rename(from = new.x) %>% 
                           rename(to = new.y) %>%
                           left_join(clade_index, by = c("from" = "node" )) %>% 
                           left_join(clade_index, by = c("to" = "node" )) %>%
                           rename(clade_from = clade.x) %>% 
                           rename(clade_to = clade.y)

# --- Compute patristic distances -------------------------------------------
pat_dist_mat <- cophenetic.phylo(tr)

# Optional: sanity checks for missing labels
missing_from <- setdiff(df_parsed$from, rownames(pat_dist_mat))
missing_to   <- setdiff(df_parsed$to,   rownames(pat_dist_mat))
if (length(missing_from) || length(missing_to)) {
  message("Warning: some labels not found in tree tip labels.\n",
          "Missing from: ", paste(missing_from, collapse = ", "),
          "\nMissing to: ", paste(missing_to, collapse = ", "))
}

# Add patristic distance for each (from, to) pair
df_parsed <- df_parsed %>% mutate(patristic = mapply(function(x, y) pat_dist_mat[x, y], from, to))
df_parsed <- df_parsed %>% mutate(id = id*100)

df_parsed <- df_parsed %>% mutate(shape = 16)
df_parsed <- df_parsed %>% mutate(size = 1.8)
df_parsed <- df_parsed %>% mutate(alpha = 0.9)
df_parsed <- df_parsed %>% mutate(color = "#696969")

#ffd454
#696969
#e62a2a

# wrong annotations
df_parsed <- df_parsed %>% mutate(color = case_when(from_superfamily != to_superfamily & clade_from == clade_to ~ "#ffd454",
                                                    from_superfamily == to_superfamily & clade_from != clade_to ~ "#ffd454",
                                                     TRUE ~ color )) %>% 
                           mutate(shape = case_when(from_superfamily != to_superfamily & clade_from == clade_to ~ 18,
                                                    from_superfamily == to_superfamily & clade_from != clade_to ~ 18,
                                                     TRUE ~ shape )) %>% 
                           mutate(size  = case_when(from_superfamily != to_superfamily & clade_from == clade_to ~ 6,
                                                    from_superfamily == to_superfamily & clade_from != clade_to ~ 6,
                                                     TRUE ~ size ))  %>% 
                           mutate(alpha  = case_when(from_superfamily != to_superfamily & clade_from == clade_to ~ 0.9,
                                                    from_superfamily == to_superfamily & clade_from != clade_to ~ 0.9,
                                                     TRUE ~ alpha ))

# clade cross events
df_parsed <- df_parsed %>% mutate(color = case_when(from_superfamily != to_superfamily & clade_from != clade_to ~ "#e62a2a",
                                                   TRUE ~ color )) %>%
                           mutate(size  = case_when(from_superfamily != to_superfamily & clade_from != clade_to ~ 12,
                                                    TRUE ~ size )) %>%
                           mutate(alpha  = case_when(from_superfamily != to_superfamily & clade_from != clade_to ~ 0.9,
                                                    TRUE ~ alpha ))

# ---------------------------------------------

first_cross <- df_parsed %>%
filter(clade_from != clade_to & from_superfamily != to_superfamily) %>%
arrange(desc(id)) %>%    # define “first” here
slice(1)

df_parsed %>%
  filter(clade_from != clade_to & from_superfamily != to_superfamily) %>%
  filter(
    (str_detect(clade_from, bad_pattern) &
      str_detect(clade_to,   bad_pattern))
  )

# --- Plot -------------------------------------------

p <- ggplot(df_parsed, aes(x = id, y = patristic, color = color)) +
      geom_point(aes(shape = shape, size = size, alpha = alpha)) +   # map inside aes()
      scale_color_identity() +
      scale_shape_identity() +
      scale_size_identity() +
      scale_alpha_identity() +
      scale_x_reverse() +
      # theme_minimal() +
      theme_linedraw() +      
      labs(x = "Percent similarity (%)", y = "Patristic Distance",
           title = "") +
      annotate("segment", x = first_cross$id, xend = first_cross$id,
                          y = 0, yend = first_cross$patristic - 0.1,
                          linetype = "dashed", color = "#484343") +
       theme(axis.title = element_text(size = 18, , margin = margin(t = 20, r = 20, b = 20, l = 20)),   # axis title size
             axis.text  = element_text(size = 16))
p
