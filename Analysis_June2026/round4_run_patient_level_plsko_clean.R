local_lib <- file.path(getwd(), "libs", paste0("R-", getRversion()[, 1], ".", getRversion()[, 2]))
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages(library(PLSKO))
suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

source("0.utils.R")

plsAKO_serial_cores1 <- function(X, y, n_ko = 25, q = 0.05, offset = 0,
                                 w.method = "lasso.lcd", covariates = NULL,
                                 gamma = 0.3, seed = 1, ...) {
  X.names <- colnames(X)
  X <- as.matrix(X)
  y <- as.vector(y)
  stopifnot(length(y) == nrow(X))
  set.seed(seed)
  pvals <- matrix(0, ncol(X), n_ko)
  selected <- vector("list", n_ko)
  for (i in seq_len(n_ko)) {
    set.seed(seed + i - 1)
    ko <- PLSKO:::plsko(X, seed = seed + i - 1, ...)
    S <- PLSKO:::ko_filter(
      X = X,
      Xk = ko,
      y = y,
      covariates = covariates,
      q = q,
      w.method = w.method,
      cores = 1
    )
    pvals[, i] <- PLSKO:::empirical_pval(S$statistic, offset = offset)
    selected[[i]] <- S
  }
  aggregated_pval <- apply(pvals, 1, PLSKO:::quantile_aggregation, gamma = gamma)
  threshold <- PLSKO:::bhq_threshold(aggregated_pval, fdr = q)
  ako.s <- which(aggregated_pval <= threshold)
  if (!is.null(X.names) && length(ako.s) > 0) {
    names(ako.s) <- X.names[ako.s]
  }
  structure(
    list(call = match.call(), s = selected, ako.s = ako.s, threshold = threshold),
    class = "AKO.result"
  )
}

input_path <- "output/plsko_patient_level/patient_level_plsko_inputs.rds"
require_analysis_files(c(
  input_path,
  "output/score_validation/round4_selected_taxa_direction_table.csv"
))

inputs_obj <- readRDS(input_path)
inputs <- inputs_obj$inputs
feature_labels <- inputs_obj$feature_labels

patient_selected <- list()
patient_status <- list()

for (label in names(inputs)) {
  message("Running clean patient-level PLSKO: ", label)
  x <- inputs[[label]]$x
  y <- inputs[[label]]$y
  
  status_row <- data.frame(
    analysis_set = label,
    n_patients = nrow(x),
    n_features = ncol(x),
    status = "not_run",
    error = NA_character_,
    n_selected = 0L
  )
  
  cache_path <- file.path("output/plsko_patient_level", paste0("clean_mean_score_", label, ".qs"))
  pls_res_ako <- tryCatch({
    if (!result_exists(cache_path)) {
      fit <- plsAKO_serial_cores1(X = x, y = y, seed = 7639)
      save_result(fit, cache_path)
      fit
    } else {
      read_result(cache_path)
    }
  }, error = function(e) {
    status_row$status <<- "failed"
    status_row$error <<- conditionMessage(e)
    NULL
  })
  
  if (is.null(pls_res_ako)) {
    patient_status[[label]] <- status_row
    next
  }
  
  selected_otu <- names(pls_res_ako$ako.s)
  status_row$status <- ifelse(length(selected_otu) > 0, "completed_selected", "completed_no_selection")
  status_row$n_selected <- length(selected_otu)
  patient_status[[label]] <- status_row
  
  if (length(selected_otu) == 0) {
    next
  }
  
  pls_res <- PhiSpace::mvr(x, y, ncomp = 2, method = "PLS", center = TRUE)
  coef_vec <- drop(pls_res$coefficients[, , 2])
  selected_labels <- feature_labels$taxon[match(selected_otu, feature_labels$otu)]
  patient_selected[[label]] <- data.frame(
    analysis_set = label,
    n_patients = nrow(x),
    n_features = ncol(x),
    otu = selected_otu,
    taxon = selected_labels,
    direction = ifelse(coef_vec[selected_otu] > 0, "positive", "negative")
  )
}

patient_selected <- bind_rows(patient_selected)
if (ncol(patient_selected) == 0) {
  patient_selected <- data.frame(
    analysis_set = character(),
    n_patients = integer(),
    n_features = integer(),
    otu = character(),
    taxon = character(),
    direction = character()
  )
}
patient_status <- bind_rows(patient_status)

write.csv(
  patient_status,
  "output/score_validation/round4_patient_level_plsko_status.csv",
  row.names = FALSE
)
write.csv(
  patient_selected,
  "output/score_validation/round4_patient_level_plsko_selected_taxa.csv",
  row.names = FALSE
)

original_direction <- read.csv(
  "output/score_validation/round4_selected_taxa_direction_table.csv",
  check.names = FALSE
)
original_long <- original_direction %>%
  pivot_longer(
    cols = c("All", "C1", "C2", "C3"),
    names_to = "analysis_set",
    values_to = "original_direction"
  ) %>%
  filter(!is.na(original_direction)) %>%
  mutate(original_direction = ifelse(original_direction == 1, "positive", "negative"))

comparison <- full_join(
  original_long,
  patient_selected %>% select(analysis_set, taxon, patient_direction = direction),
  by = c("analysis_set", "taxon")
) %>%
  mutate(
    selected_original = !is.na(original_direction),
    selected_patient_level = !is.na(patient_direction),
    direction_concordant = selected_original & selected_patient_level &
      original_direction == patient_direction
  )

write.csv(
  comparison,
  "output/score_validation/round4_patient_vs_timepoint_plsko_comparison.csv",
  row.names = FALSE
)

comparison_summary <- comparison %>%
  group_by(analysis_set) %>%
  summarise(
    n_original = sum(selected_original),
    n_patient_level = sum(selected_patient_level),
    n_overlap = sum(selected_original & selected_patient_level),
    n_direction_concordant = sum(direction_concordant),
    original_taxa = paste(taxon[selected_original], collapse = "; "),
    patient_level_taxa = paste(taxon[selected_patient_level], collapse = "; "),
    overlap_taxa = paste(taxon[selected_original & selected_patient_level], collapse = "; "),
    .groups = "drop"
  )

write.csv(
  comparison_summary,
  "output/score_validation/round4_patient_vs_timepoint_plsko_summary.csv",
  row.names = FALSE
)

print(patient_status)
print(comparison_summary)
