local_lib <- file.path(getwd(), "libs", paste0("R-", getRversion()[, 1], ".", getRversion()[, 2]))
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(nlme))

source("0.utils.R")

dir.create("output/score_validation", recursive = TRUE, showWarnings = FALSE)

require_analysis_files(c(
  "data/patientdata_3obs_allresponses.csv",
  "data/OTU_table.csv",
  "data/taxonomy_table.csv",
  "data_filtered/meta_filtered.csv",
  "output/score_validation/round4_cluster_assignments_from_render.csv",
  "output/score_validation/round4_selected_taxa_direction_table.csv",
  "output/score_validation/round4_repeated_measure_lmm_taxa_screen.csv"
))

fit_lmm_screen <- function(clr_sce, screen_meta, taxon_labels) {
  sample_ids <- intersect(rownames(screen_meta), colnames(clr_sce))
  screen_meta <- screen_meta[sample_ids, c("id", "timepoint", "cluster", "OMscore_calib")]
  screen_x <- t(assay(clr_sce, "data"))[sample_ids, , drop = FALSE]

  bind_rows(lapply(colnames(screen_x), function(otu) {
    model_df <- data.frame(
      abundance = as.numeric(screen_x[, otu]),
      screen_meta
    )
    fit <- tryCatch(
      nlme::lme(
        abundance ~ OMscore_calib + timepoint + cluster,
        random = ~1 | id,
        data = model_df,
        na.action = na.omit,
        control = nlme::lmeControl(opt = "optim")
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      return(data.frame(
        otu = otu,
        taxon = taxon_labels[[otu]],
        beta_om_score = NA_real_,
        p_om_score = NA_real_,
        error = conditionMessage(fit)
      ))
    }
    tab <- summary(fit)$tTable
    data.frame(
      otu = otu,
      taxon = taxon_labels[[otu]],
      beta_om_score = tab["OMscore_calib", "Value"],
      p_om_score = tab["OMscore_calib", "p-value"],
      error = NA_character_
    )
  })) %>%
    mutate(
      fdr_om_score = p.adjust(p_om_score, method = "BH"),
      direction = ifelse(beta_om_score > 0, "positive", "negative")
    ) %>%
    arrange(fdr_om_score, p_om_score)
}

meta_dat <- read.csv("data_filtered/meta_filtered.csv", row.names = 1) %>%
  select(c(id, distress:teeth, timepoint))
symptom_matrix <- meta_dat %>% select(distress:teeth) %>% as.matrix()
for (icol in which(colSums(is.na(symptom_matrix)) > 0)) {
  vals <- symptom_matrix[, icol]
  symptom_matrix[is.na(vals), icol] <- mean(vals, na.rm = TRUE)
}
set.seed(8790)
npca <- nsprcomp::nsprcomp(PhiSpace::doubleCent(symptom_matrix), center = FALSE, nneg = TRUE)
meta_dat$OMscore_calib <- npca$x[, "PC1"]

cluster_assignments <- read.csv("output/score_validation/round4_cluster_assignments_from_render.csv")
meta_dat$cluster <- factor(cluster_assignments$cluster[match(meta_dat$id, cluster_assignments$id)])

taxonomy_table <- read.csv("data/taxonomy_table.csv", row.names = 1)
taxon_genus <- taxonomy_table[, "genus"]
taxon_genus[is.na(taxon_genus) | taxon_genus == ""] <- "Unclassified"
taxon_labels <- setNames(
  paste0(taxon_genus, gsub("OTU", "", rownames(taxonomy_table))),
  rownames(taxonomy_table)
)

otu_raw <- read.csv("data/OTU_table.csv", row.names = 1, check.names = FALSE)
new_names <- sub(".*MB", "MB", colnames(otu_raw))
new_names <- sub(".*MF", "MF", new_names)
colnames(otu_raw) <- new_names
common_samples <- intersect(new_names, rownames(meta_dat))
otu_counts <- t(otu_raw[, common_samples, drop = FALSE])

thresholds <- c(0.005, 0.01, 0.02)
current_threshold <- 0.01

original_direction <- read.csv(
  "output/score_validation/round4_selected_taxa_direction_table.csv",
  check.names = FALSE
)
original_selected_taxa <- original_direction %>%
  filter(if_any(-taxon, ~ !is.na(.x))) %>%
  pull(taxon) %>%
  unique()
selected_otu_by_taxon <- names(taxon_labels)[match(original_selected_taxa, taxon_labels)]
selected_otu_by_taxon <- selected_otu_by_taxon[!is.na(selected_otu_by_taxon)]

reference_lmm <- read.csv("output/score_validation/round4_repeated_measure_lmm_taxa_screen.csv") %>%
  filter(fdr_om_score < 0.05)
reference_lmm_otus <- reference_lmm$otu

threshold_results <- list()
screen_results <- list()

for (threshold in thresholds) {
  rel_abund_percent <- sweep(otu_counts, 1, rowSums(otu_counts), FUN = "/") * 100
  mean_ra <- colMeans(rel_abund_percent, na.rm = TRUE)
  keep_otus <- names(mean_ra)[mean_ra >= threshold]
  filtered_counts <- otu_counts[, keep_otus, drop = FALSE]
  filtered_sce <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = t(filtered_counts)),
    colData = meta_dat[rownames(filtered_counts), , drop = FALSE]
  )
  filtered_clr <- PhiSpace::CLRnorm(filtered_sce)
  lmm_results <- fit_lmm_screen(filtered_clr, meta_dat, taxon_labels)
  threshold_tag <- gsub("\\.", "p", sprintf("%.3f", threshold))

  write.csv(
    lmm_results,
    file.path(
      "output/score_validation",
      paste0("round4_filter_", threshold_tag, "_lmm_taxa_screen.csv")
    ),
    row.names = FALSE
  )

  screen_results[[as.character(threshold)]] <- lmm_results
  threshold_results[[as.character(threshold)]] <- data.frame(
    mean_ra_threshold_percent = threshold,
    n_retained_otus = length(keep_otus),
    total_count_mass_removed_percent = (1 - sum(filtered_counts, na.rm = TRUE) / sum(otu_counts, na.rm = TRUE)) * 100,
    n_current_147_retained = sum(names(mean_ra)[mean_ra >= current_threshold] %in% keep_otus),
    n_original_plsko_selected_retained = sum(selected_otu_by_taxon %in% keep_otus),
    n_original_plsko_selected_total = length(selected_otu_by_taxon),
    n_lmm_fdr_0_05 = sum(lmm_results$fdr_om_score < 0.05, na.rm = TRUE),
    n_reference_lmm_fdr_0_05_retained = sum(reference_lmm_otus %in% keep_otus),
    n_reference_lmm_fdr_0_05_still_significant = sum(
      lmm_results$otu %in% reference_lmm_otus & lmm_results$fdr_om_score < 0.05,
      na.rm = TRUE
    )
  )
}

summary_df <- bind_rows(threshold_results)
write.csv(
  summary_df,
  "output/score_validation/round4_filter_threshold_sensitivity_summary.csv",
  row.names = FALSE
)

reference_threshold_results <- screen_results[[as.character(current_threshold)]] %>%
  select(
    otu,
    taxon,
    current_beta = beta_om_score,
    current_fdr = fdr_om_score,
    current_direction = direction
  )

comparison_df <- bind_rows(lapply(names(screen_results), function(threshold) {
  screen_results[[threshold]] %>%
    mutate(mean_ra_threshold_percent = as.numeric(threshold)) %>%
    select(
      mean_ra_threshold_percent,
      otu,
      taxon,
      beta_om_score,
      fdr_om_score,
      direction
    )
})) %>%
  left_join(reference_threshold_results, by = c("otu", "taxon")) %>%
  mutate(
    reference_lmm_fdr_0_05 = !is.na(current_fdr) & current_fdr < 0.05,
    threshold_lmm_fdr_0_05 = !is.na(fdr_om_score) & fdr_om_score < 0.05,
    direction_concordant_with_current = direction == current_direction
  ) %>%
  arrange(mean_ra_threshold_percent, fdr_om_score, otu)

write.csv(
  comparison_df,
  "output/score_validation/round4_filter_threshold_lmm_comparison.csv",
  row.names = FALSE
)

selected_retention <- expand.grid(
  mean_ra_threshold_percent = thresholds,
  otu = selected_otu_by_taxon,
  stringsAsFactors = FALSE
) %>%
  mutate(taxon = taxon_labels[otu]) %>%
  left_join(
    bind_rows(lapply(thresholds, function(threshold) {
      rel_abund_percent <- sweep(otu_counts, 1, rowSums(otu_counts), FUN = "/") * 100
      mean_ra <- colMeans(rel_abund_percent, na.rm = TRUE)
      data.frame(
        mean_ra_threshold_percent = threshold,
        otu = names(mean_ra),
        mean_relative_abundance_percent = mean_ra,
        retained = mean_ra >= threshold
      )
    })),
    by = c("mean_ra_threshold_percent", "otu")
  ) %>%
  left_join(
    comparison_df %>%
      select(
        mean_ra_threshold_percent,
        otu,
        beta_om_score,
        fdr_om_score,
        direction,
        threshold_lmm_fdr_0_05
      ),
    by = c("mean_ra_threshold_percent", "otu")
  ) %>%
  arrange(mean_ra_threshold_percent, taxon)

write.csv(
  selected_retention,
  "output/score_validation/round4_filter_threshold_selected_taxa_retention.csv",
  row.names = FALSE
)

print(summary_df)
