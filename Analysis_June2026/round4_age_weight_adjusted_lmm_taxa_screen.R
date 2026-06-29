local_lib <- file.path(getwd(), "libs", paste0("R-", getRversion()[, 1], ".", getRversion()[, 2]))
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(dplyr))
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
  "output/score_validation/round4_selected_taxa_direction_table.csv",
  "output/score_validation/round4_repeated_measure_lmm_taxa_screen.csv",
  "output/score_validation/round4_demographic_merged_readr.csv"
))

to_numeric <- function(x) {
  as.numeric(gsub("[^0-9.\\-]", "", as.character(x)))
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
meta_dat$cluster <- cluster_assignments$cluster[match(meta_dat$id, cluster_assignments$id)]
meta_dat$cluster <- factor(meta_dat$cluster)

demographics <- read.csv("output/score_validation/round4_demographic_merged_readr.csv") %>%
  transmute(
    id,
    age = to_numeric(age),
    weight = to_numeric(weight)
  )

meta_dat <- meta_dat %>%
  left_join(demographics, by = "id") %>%
  mutate(
    age_z = as.numeric(scale(age)),
    weight_z = as.numeric(scale(weight))
  )
rownames(meta_dat) <- rownames(read.csv("data_filtered/meta_filtered.csv", row.names = 1))

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
screen_meta <- meta_dat[sample_ids, c("id", "timepoint", "cluster", "OMscore_calib", "age_z", "weight_z")]
complete_sample_ids <- rownames(screen_meta)[complete.cases(screen_meta)]
screen_meta <- screen_meta[complete_sample_ids, , drop = FALSE]
screen_x <- t(assay(microb, "data"))[complete_sample_ids, , drop = FALSE]

lmm_results <- bind_rows(lapply(colnames(screen_x), function(otu) {
  model_df <- data.frame(
    abundance = as.numeric(screen_x[, otu]),
    screen_meta
  )
  fit <- tryCatch(
    nlme::lme(
      abundance ~ OMscore_calib + timepoint + cluster + age_z + weight_z,
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
      beta_om_score_adjusted = NA_real_,
      p_om_score_adjusted = NA_real_,
      beta_age_z = NA_real_,
      p_age_z = NA_real_,
      beta_weight_z = NA_real_,
      p_weight_z = NA_real_,
      error = conditionMessage(fit)
    ))
  }
  tab <- summary(fit)$tTable
  data.frame(
    otu = otu,
    taxon = taxon_labels[[otu]],
    beta_om_score_adjusted = tab["OMscore_calib", "Value"],
    p_om_score_adjusted = tab["OMscore_calib", "p-value"],
    beta_age_z = tab["age_z", "Value"],
    p_age_z = tab["age_z", "p-value"],
    beta_weight_z = tab["weight_z", "Value"],
    p_weight_z = tab["weight_z", "p-value"],
    error = NA_character_
  )
}))

lmm_results <- lmm_results %>%
  mutate(
    fdr_om_score_adjusted = p.adjust(p_om_score_adjusted, method = "BH"),
    fdr_age_z = p.adjust(p_age_z, method = "BH"),
    fdr_weight_z = p.adjust(p_weight_z, method = "BH"),
    direction_adjusted = ifelse(beta_om_score_adjusted > 0, "positive", "negative")
  ) %>%
  arrange(fdr_om_score_adjusted, p_om_score_adjusted)

write.csv(
  lmm_results,
  "output/score_validation/round4_age_weight_adjusted_lmm_taxa_screen.csv",
  row.names = FALSE
)

unadjusted <- read.csv("output/score_validation/round4_repeated_measure_lmm_taxa_screen.csv") %>%
  select(
    otu,
    taxon,
    beta_om_score_unadjusted = beta_om_score,
    p_om_score_unadjusted = p_om_score,
    fdr_om_score_unadjusted = fdr_om_score,
    direction_unadjusted = direction
  )

comparison <- unadjusted %>%
  full_join(
    lmm_results %>%
      select(
        otu,
        taxon,
        beta_om_score_adjusted,
        p_om_score_adjusted,
        fdr_om_score_adjusted,
        direction_adjusted
      ),
    by = c("otu", "taxon")
  ) %>%
  mutate(
    unadjusted_fdr_0_05 = !is.na(fdr_om_score_unadjusted) & fdr_om_score_unadjusted < 0.05,
    adjusted_fdr_0_05 = !is.na(fdr_om_score_adjusted) & fdr_om_score_adjusted < 0.05,
    direction_concordant = direction_unadjusted == direction_adjusted
  ) %>%
  arrange(fdr_om_score_adjusted, p_om_score_adjusted)

write.csv(
  comparison,
  "output/score_validation/round4_age_weight_adjusted_lmm_comparison.csv",
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

adjusted_original_overlap <- original_global %>%
  left_join(lmm_results, by = "taxon") %>%
  mutate(
    adjusted_lmm_fdr_0_05 = !is.na(fdr_om_score_adjusted) & fdr_om_score_adjusted < 0.05,
    direction_concordant = original_direction == direction_adjusted
  ) %>%
  arrange(fdr_om_score_adjusted, p_om_score_adjusted)

write.csv(
  adjusted_original_overlap,
  "output/score_validation/round4_age_weight_adjusted_lmm_original_plsko_overlap.csv",
  row.names = FALSE
)

summary_df <- data.frame(
  n_samples_complete = nrow(screen_meta),
  n_patients_complete = length(unique(screen_meta$id)),
  n_taxa_tested = nrow(lmm_results),
  n_adjusted_lmm_fdr_0_05 = sum(lmm_results$fdr_om_score_adjusted < 0.05, na.rm = TRUE),
  n_unadjusted_lmm_fdr_0_05 = sum(comparison$unadjusted_fdr_0_05, na.rm = TRUE),
  n_unadjusted_supported_still_adjusted_fdr_0_05 = sum(
    comparison$unadjusted_fdr_0_05 & comparison$adjusted_fdr_0_05,
    na.rm = TRUE
  ),
  n_unadjusted_supported_direction_concordant = sum(
    comparison$unadjusted_fdr_0_05 & comparison$direction_concordant,
    na.rm = TRUE
  ),
  n_original_global_plsko = nrow(original_global),
  n_original_global_adjusted_lmm_fdr_0_05 = sum(
    adjusted_original_overlap$adjusted_lmm_fdr_0_05,
    na.rm = TRUE
  ),
  n_original_global_adjusted_direction_concordant = sum(
    adjusted_original_overlap$adjusted_lmm_fdr_0_05 &
      adjusted_original_overlap$direction_concordant,
    na.rm = TRUE
  )
)

write.csv(
  summary_df,
  "output/score_validation/round4_age_weight_adjusted_lmm_summary.csv",
  row.names = FALSE
)

print(summary_df)
print(head(lmm_results, 20))
