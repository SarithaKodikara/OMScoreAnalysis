local_lib <- file.path(getwd(), "libs", paste0("R-", getRversion()[, 1], ".", getRversion()[, 2]))
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(SummarizedExperiment))

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
  "output/score_validation/round4_age_weight_adjusted_lmm_taxa_screen.csv"
))

if (!result_exists("output/plsko/global.qs")) {
  stop("Missing cached global PLSKO result.", call. = FALSE)
}
if (!result_exists("output/plsko/group_list.qs")) {
  stop("Missing cached cluster-specific PLSKO result.", call. = FALSE)
}

direction_label <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x == 1 ~ "positive",
    x == 0 ~ "negative",
    TRUE ~ NA_character_
  )
}

clean_taxonomy <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unclassified"
  x
}

get_pls_coefficients <- function(x, y, taxon_labels) {
  keep <- colVars(x, na.rm = TRUE) > 1e-10
  x <- x[, keep, drop = FALSE]
  fit <- PhiSpace::mvr(x, y, ncomp = 2, method = "PLS", center = TRUE)
  coef_vec <- as.numeric(fit$coefficients[, , 2])
  names(coef_vec) <- rownames(fit$coefficients)
  tibble(
    otu = names(coef_vec),
    taxon = unname(taxon_labels[names(coef_vec)]),
    pls_coefficient = coef_vec
  )
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

taxonomy_table <- read.csv("data/taxonomy_table.csv", row.names = 1)
taxonomy_table <- taxonomy_table %>%
  mutate(across(everything(), clean_taxonomy)) %>%
  tibble::rownames_to_column("otu")

taxon_labels <- setNames(
  paste0(taxonomy_table$genus, gsub("OTU", "", taxonomy_table$otu)),
  taxonomy_table$otu
)

microb <- readRDS("data_filtered/sce_all.rds")
selected_otu <- readRDS("data_filtered/selectedOTUs.rds")
microb <- PhiSpace::CLRnorm(microb[selected_otu, ])

sample_ids <- intersect(rownames(meta_dat), colnames(microb))
meta_dat <- meta_dat[sample_ids, , drop = FALSE]
microb_mat <- t(assay(microb, "data"))[sample_ids, , drop = FALSE]

global_plsko <- read_result("output/plsko/global.qs")
group_plsko <- read_result("output/plsko/group_list.qs")

direction_table <- read.csv(
  "output/score_validation/round4_selected_taxa_direction_table.csv",
  check.names = FALSE
)

selected_labels <- unique(direction_table$taxon)
label_to_otu <- tibble(
  otu = names(taxon_labels),
  taxon = unname(taxon_labels)
) %>%
  filter(taxon %in% selected_labels)

global_coef <- get_pls_coefficients(
  microb_mat,
  meta_dat$OMscore_calib,
  taxon_labels
) %>%
  transmute(taxon, pls_coefficient_global = pls_coefficient)

cluster_coef <- bind_rows(lapply(seq_along(group_plsko), function(cluster_id) {
  patient_ids <- unique(meta_dat$id[meta_dat$cluster == cluster_id])
  cluster_samples <- rownames(meta_dat)[meta_dat$id %in% patient_ids]
  get_pls_coefficients(
    microb_mat[cluster_samples, , drop = FALSE],
    meta_dat[cluster_samples, "OMscore_calib"],
    taxon_labels
  ) %>%
    mutate(cluster = paste0("C", cluster_id))
})) %>%
  select(taxon, cluster, pls_coefficient) %>%
  pivot_wider(
    names_from = cluster,
    values_from = pls_coefficient,
    names_prefix = "pls_coefficient_"
  )

lmm_primary <- read.csv("output/score_validation/round4_repeated_measure_lmm_taxa_screen.csv") %>%
  select(
    otu,
    lmm_beta_om_score = beta_om_score,
    lmm_p_om_score = p_om_score,
    lmm_fdr_om_score = fdr_om_score,
    lmm_direction = direction
  )

lmm_adjusted <- read.csv("output/score_validation/round4_age_weight_adjusted_lmm_taxa_screen.csv") %>%
  select(
    otu,
    lmm_adjusted_beta_om_score = beta_om_score_adjusted,
    lmm_adjusted_p_om_score = p_om_score_adjusted,
    lmm_adjusted_fdr_om_score = fdr_om_score_adjusted,
    lmm_adjusted_direction = direction_adjusted
  )

supplementary_table <- direction_table %>%
  mutate(
    direction_global = direction_label(All),
    direction_C1 = direction_label(C1),
    direction_C2 = direction_label(C2),
    direction_C3 = direction_label(C3),
    selected_global = !is.na(All),
    selected_C1 = !is.na(C1),
    selected_C2 = !is.na(C2),
    selected_C3 = !is.na(C3)
  ) %>%
  select(
    taxon,
    selected_global, direction_global,
    selected_C1, direction_C1,
    selected_C2, direction_C2,
    selected_C3, direction_C3
  ) %>%
  left_join(label_to_otu, by = "taxon") %>%
  left_join(taxonomy_table, by = "otu") %>%
  left_join(global_coef, by = "taxon") %>%
  left_join(cluster_coef, by = "taxon") %>%
  left_join(lmm_primary, by = "otu") %>%
  left_join(lmm_adjusted, by = "otu") %>%
  mutate(
    pls_coefficient_global = ifelse(
      selected_global,
      pls_coefficient_global,
      NA_real_
    ),
    pls_coefficient_C1 = ifelse(selected_C1, pls_coefficient_C1, NA_real_),
    pls_coefficient_C2 = ifelse(selected_C2, pls_coefficient_C2, NA_real_),
    pls_coefficient_C3 = ifelse(selected_C3, pls_coefficient_C3, NA_real_),
    lmm_fdr_lt_0_05 = !is.na(lmm_fdr_om_score) & lmm_fdr_om_score < 0.05,
    lmm_adjusted_fdr_lt_0_05 = !is.na(lmm_adjusted_fdr_om_score) &
      lmm_adjusted_fdr_om_score < 0.05
  ) %>%
  select(
    otu,
    taxon,
    kingdom,
    phylum,
    class,
    order,
    family,
    genus,
    selected_global,
    direction_global,
    pls_coefficient_global,
    selected_C1,
    direction_C1,
    pls_coefficient_C1,
    selected_C2,
    direction_C2,
    pls_coefficient_C2,
    selected_C3,
    direction_C3,
    pls_coefficient_C3,
    lmm_beta_om_score,
    lmm_p_om_score,
    lmm_fdr_om_score,
    lmm_direction,
    lmm_fdr_lt_0_05,
    lmm_adjusted_beta_om_score,
    lmm_adjusted_p_om_score,
    lmm_adjusted_fdr_om_score,
    lmm_adjusted_direction,
    lmm_adjusted_fdr_lt_0_05
  ) %>%
  arrange(genus, otu)

write.csv(
  supplementary_table,
  "output/score_validation/round4_selected_taxa_supplementary_table.csv",
  row.names = FALSE
)

summary_table <- data.frame(
  n_selected_taxa = nrow(supplementary_table),
  n_global_selected = sum(supplementary_table$selected_global),
  n_cluster1_selected = sum(supplementary_table$selected_C1),
  n_cluster2_selected = sum(supplementary_table$selected_C2),
  n_cluster3_selected = sum(supplementary_table$selected_C3),
  n_primary_lmm_fdr_0_05 = sum(supplementary_table$lmm_fdr_lt_0_05, na.rm = TRUE),
  n_adjusted_lmm_fdr_0_05 = sum(supplementary_table$lmm_adjusted_fdr_lt_0_05, na.rm = TRUE)
)

write.csv(
  summary_table,
  "output/score_validation/round4_selected_taxa_supplementary_table_summary.csv",
  row.names = FALSE
)

print(summary_table)
print(supplementary_table)
