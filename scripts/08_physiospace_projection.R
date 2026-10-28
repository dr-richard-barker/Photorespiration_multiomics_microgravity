#!/usr/bin/env Rscript
# 08_physiospace_projection.R
# Calculates PhysioSpace stress pattern scores for all studies in the registry.
#
# Requirements:
#   remotes::install_git("https://git.rwth-aachen.de/jrc-combine/PhysioSpaceMethods")
#   remotes::install_git("https://git.rwth-aachen.de/jrc-combine/PlantPhysioSpace")
#
# Output:
#   results/tables/T14_physiospace_scores.csv
#
# STATUS (2026-10-27): this script has never run to completion. The packages are not
# installed here, and `calculatePhysioScores` / `AT_Stress_Space` have not been checked
# against the packages' real exports. Any T14 file that predates a real run is not output
# of this script.

library(dplyr)
library(tidyr)
library(readr)
library(PhysioSpaceMethods)
library(PlantPhysioSpace)

cat("Loading registry...\n")
registry <- read_tsv("data/study_registry.tsv", comment = "#", show_col_types = FALSE)
included <- registry %>% filter(include == "yes")

# We will collect Entrez-mapped Log2FC from all cache files
# For simplicity, assuming data/cache/dge_{accession}.csv format,
# and that ENTREZID is present.

all_data <- list()
for (i in 1:nrow(included)) {
    acc <- included$accession[i]
    contrast_name <- included$contrast[i]
    base_acc <- strsplit(strsplit(acc, "-light")[[1]][1], "-dark")[[1]][1]
    
    file_path <- file.path("data", "cache", paste0("dge_", base_acc, ".csv"))
    if (!file.exists(file_path)) {
        file_path <- file.path("data", "cache", included$dge_file[i])
    }
    
    if (file.exists(file_path)) {
        df <- read_csv(file_path, show_col_types = FALSE)
        if ("ENTREZID" %in% names(df) && contrast_name %in% names(df)) {
            df_sub <- df %>% 
                select(ENTREZID, all_of(contrast_name)) %>%
                rename(!!acc := all_of(contrast_name)) %>%
                filter(!is.na(ENTREZID), !is.na(get(acc))) %>%
                mutate(ENTREZID = as.character(ENTREZID)) %>%
                group_by(ENTREZID) %>% summarize_all(mean) # handle duplicates
            all_data[[acc]] <- df_sub
        }
    }
}

missing <- setdiff(included$accession, names(all_data))
if (length(missing) > 0) {
    stop("No usable log2FC (ENTREZID + contrast column) for: ", paste(missing, collapse = ", "))
}

cat("Merging data...\n")
combined_df <- Reduce(function(x, y) full_join(x, y, by = "ENTREZID"), all_data)
combined_df[is.na(combined_df)] <- 0

# Convert to matrix
expr_mat <- as.matrix(combined_df[,-1])
rownames(expr_mat) <- combined_df$ENTREZID

cat("Running Plant PhysioSpace...\n")
# Load the Arabidopsis stress space from the package
data("AT_Stress_Space", package="PlantPhysioSpace")

# Calculate PhysioScores
# Default method uses dot product or correlation. PlantPhysioSpace often uses calculatePhysioScores.
scores <- calculatePhysioScores(expr_mat, Space = AT_Stress_Space)

cat("Saving results...\n")
out_dir <- "results/tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(scores, file.path(out_dir, "T14_physiospace_scores.csv"), row.names = TRUE)

cat("Done! Output saved to T14_physiospace_scores.csv\n")
