local_lib <- file.path(getwd(), "libs", paste0("R-", getRversion()[, 1], ".", getRversion()[, 2]))
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(nlme))

source("0.utils.R")

dir.create("output/score_validation", recursive = TRUE, showWarnings = FALSE)

require_analysis_files(c(
  "data_filtered/meta_filtered.csv",
  "data/taxonomy_table.csv",
  "data_filtered/sce_all.rds",
  "data_filtered/selectedOTUs.rds",
  "output/score_validation/round4_cluster_assignments_from_render.csv",
  "output/score_validation/round4_selected_taxa_direction_table.csv"
))

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
meta_dat$cluster <- cluster_assignments$cluster[match(meta_dat$id, cluster_assignments$id)]
meta_dat$cluster <- factor(meta_dat$cluster)

taxonomy_table <- read.csv("data/taxonomy_table.csv", row.names = 1)
taxon_genus <- taxonomy_table[, "genus"]
taxon_genus[is.na(taxon_genus) | taxon_genus == ""] <- "Unclassified"
taxon_labels <- setNames(
  paste0(taxon_genus, gsub("OTU", "", rownames(taxonomy_table))),
  rownames(taxonomy_table)
)

microb <- readRDS("data_filtered/sce_all.rds")
selected_otu <- readRDS("data_filtered/selectedOTUs.rds")
microb <- PhiSpace::CLRnorm(microb[selected_otu, ])

sample_ids <- intersect(rownames(meta_dat), colnames(microb))
screen_meta <- meta_dat[sample_ids, c("id", "timepoint", "cluster", "OMscore_calib")]
screen_x <- t(assay(microb, "data"))[sample_ids, , drop = FALSE]
lmm_results <- bind_rows(lapply(colnames(screen_x), function(otu) {
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
}))

lmm_results <- lmm_results %>%
  mutate(
    fdr_om_score = p.adjust(p_om_score, method = "BH"),
    direction = ifelse(beta_om_score > 0, "positive", "negative")
  ) %>%
  arrange(fdr_om_score, p_om_score)

write.csv(
  lmm_results,
  "output/score_validation/round4_repeated_measure_lmm_taxa_screen.csv",
  row.names = FALSE
)

original_direction <- read.csv(
  "output/score_validation/round4_selected_taxa_direction_table.csv",
  check.names = FALSE
)
original_global <- original_direction %>%
  filter(!is.na(All)) %>%
  transmute(
    taxon,
    original_direction = ifelse(All == 1, "positive", "negative")
  )

lmm_original_overlap <- original_global %>%
  left_join(lmm_results, by = "taxon") %>%
  mutate(
    lmm_fdr_0_05 = !is.na(fdr_om_score) & fdr_om_score < 0.05,
    direction_concordant = original_direction == direction
  ) %>%
  arrange(fdr_om_score, p_om_score)

write.csv(
  lmm_original_overlap,
  "output/score_validation/round4_repeated_measure_lmm_original_plsko_overlap.csv",
  row.names = FALSE
)

lmm_summary <- data.frame(
  n_taxa_tested = nrow(lmm_results),
  n_lmm_fdr_0_05 = sum(lmm_results$fdr_om_score < 0.05, na.rm = TRUE),
  n_original_global_plsko = nrow(original_global),
  n_original_global_lmm_fdr_0_05 = sum(lmm_original_overlap$lmm_fdr_0_05, na.rm = TRUE),
  n_original_global_direction_concordant = sum(
    lmm_original_overlap$lmm_fdr_0_05 & lmm_original_overlap$direction_concordant,
    na.rm = TRUE
  )
)

write.csv(
  lmm_summary,
  "output/score_validation/round4_repeated_measure_lmm_summary.csv",
  row.names = FALSE
)

print(lmm_summary)
print(head(lmm_results, 20))
