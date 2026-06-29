local_lib <- file.path(getwd(), "libs", paste0("R-", getRversion()[, 1], ".", getRversion()[, 2]))
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(SummarizedExperiment))

source("0.utils.R")

dir.create("output/plsko_patient_level", recursive = TRUE, showWarnings = FALSE)

require_analysis_files(c(
  "data_filtered/meta_filtered.csv",
  "data/taxonomy_table.csv",
  "data_filtered/sce_all.rds",
  "data_filtered/selectedOTUs.rds",
  "output/score_validation/round4_cluster_assignments_from_render.csv"
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
taxonomy_table <- read.csv("data/taxonomy_table.csv", row.names = 1)
microb <- readRDS("data_filtered/sce_all.rds")
selected_otu <- readRDS("data_filtered/selectedOTUs.rds")
microb <- PhiSpace::CLRnorm(microb[selected_otu, ])

sample_ids <- intersect(rownames(meta_dat), colnames(microb))
sample_meta <- meta_dat[sample_ids, c("id", "timepoint", "OMscore_calib")]
sample_x <- t(assay(microb, "data"))[sample_ids, , drop = FALSE]

patient_counts <- table(sample_meta$id)
patient_x_sum <- rowsum(sample_x, group = sample_meta$id)
patient_x <- sweep(
  patient_x_sum,
  1,
  as.numeric(patient_counts[rownames(patient_x_sum)]),
  FUN = "/"
)

patient_y <- sample_meta %>%
  group_by(id) %>%
  summarise(
    mean_score = mean(OMscore_calib, na.rm = TRUE),
    max_score = max(OMscore_calib, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  left_join(cluster_assignments, by = "id")

patient_x <- patient_x[patient_y$id, , drop = FALSE]
stopifnot(identical(rownames(patient_x), patient_y$id))

analysis_sets <- list(
  All = seq_len(nrow(patient_y)),
  C1 = which(patient_y$cluster == 1),
  C2 = which(patient_y$cluster == 2),
  C3 = which(patient_y$cluster == 3)
)

inputs <- lapply(names(analysis_sets), function(label) {
  idx <- analysis_sets[[label]]
  x <- as.matrix(patient_x[idx, , drop = FALSE])
  y <- as.numeric(patient_y$mean_score[idx])
  keep <- colVars(x) > 1e-10
  list(
    label = label,
    x = x[, keep, drop = FALSE],
    y = y,
    patient_ids = rownames(x)
  )
})
names(inputs) <- names(analysis_sets)

taxon_genus <- taxonomy_table[, "genus"]
taxon_genus[is.na(taxon_genus) | taxon_genus == ""] <- "Unclassified"
feature_labels <- data.frame(
  otu = rownames(taxonomy_table),
  taxon = paste0(taxon_genus, gsub("OTU", "", rownames(taxonomy_table)))
)

saveRDS(
  list(
    inputs = inputs,
    patient_y = patient_y,
    feature_labels = feature_labels
  ),
  "output/plsko_patient_level/patient_level_plsko_inputs.rds"
)

write.csv(
  patient_y,
  "output/score_validation/round4_patient_level_plsko_patient_inputs.csv",
  row.names = FALSE
)

cat("Saved patient-level PLSKO inputs for", length(inputs), "analysis sets.\n")
