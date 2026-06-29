local_lib <- file.path(getwd(), "libs", paste0("R-", getRversion()[, 1], ".", getRversion()[, 2]))
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages(library(dplyr))

source("0.utils.R")

dir.create("output/score_validation", recursive = TRUE, showWarnings = FALSE)

require_analysis_files(c(
  "data/age_final_data.csv",
  "output/score_validation/round4_cluster_assignments_from_render.csv"
))

to_numeric <- function(x) {
  as.numeric(gsub("[^0-9.\\-]", "", as.character(x)))
}

pooled_sd <- function(x, y) {
  nx <- sum(!is.na(x))
  ny <- sum(!is.na(y))
  sqrt(((nx - 1) * var(x, na.rm = TRUE) + (ny - 1) * var(y, na.rm = TRUE)) / (nx + ny - 2))
}

cohen_d <- function(x, y) {
  (mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)) / pooled_sd(x, y)
}

cluster_assignments <- read.csv("output/score_validation/round4_cluster_assignments_from_render.csv")

demographics <- readr::read_csv("data/age_final_data.csv", show_col_types = FALSE) %>%
  right_join(cluster_assignments, by = "id") %>%
  transmute(
    id,
    cluster = factor(cluster),
    age = to_numeric(age),
    weight_lb = to_numeric(weight),
    weight_kg = weight_lb * 0.45359237,
    height_in = to_numeric(heightnew),
    height_cm = height_in * 2.54,
    sex = factor(sex, levels = c(1, 2), labels = c("Female", "Male")),
    is_non_hispanic_white = race___1 == 1
  )

summary_table <- demographics %>%
  group_by(cluster) %>%
  summarise(
    n_cluster = n(),
    n_demo = sum(complete.cases(age)),
    age_mean = mean(age, na.rm = TRUE),
    age_sd = sd(age, na.rm = TRUE),
    weight_lb_mean = mean(weight_lb, na.rm = TRUE),
    weight_lb_sd = sd(weight_lb, na.rm = TRUE),
    weight_kg_mean = mean(weight_kg, na.rm = TRUE),
    weight_kg_sd = sd(weight_kg, na.rm = TRUE),
    height_in_mean = mean(height_in, na.rm = TRUE),
    height_in_sd = sd(height_in, na.rm = TRUE),
    height_cm_mean = mean(height_cm, na.rm = TRUE),
    height_cm_sd = sd(height_cm, na.rm = TRUE),
    male = sum(sex == "Male", na.rm = TRUE),
    female = sum(sex == "Female", na.rm = TRUE),
    non_hispanic_white = sum(is_non_hispanic_white, na.rm = TRUE),
    other_ethnicity = sum(!is_non_hispanic_white & !is.na(is_non_hispanic_white), na.rm = TRUE),
    .groups = "drop"
  )

pairs <- combn(levels(demographics$cluster), 2, simplify = FALSE)
variables <- c(age = "age", weight = "weight_lb", height = "height_in")

pairwise_results <- bind_rows(lapply(names(variables), function(variable_label) {
  variable <- variables[[variable_label]]
  bind_rows(lapply(pairs, function(pair) {
    x <- demographics[[variable]][demographics$cluster == pair[1]]
    y <- demographics[[variable]][demographics$cluster == pair[2]]
    test <- t.test(x, y, var.equal = TRUE)
    data.frame(
      variable = variable_label,
      comparison = paste0("C", pair[1], " vs C", pair[2]),
      cluster_1 = pair[1],
      cluster_2 = pair[2],
      n_1 = sum(!is.na(x)),
      n_2 = sum(!is.na(y)),
      mean_1 = mean(x, na.rm = TRUE),
      mean_2 = mean(y, na.rm = TRUE),
      mean_difference = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
      cohen_d = cohen_d(x, y),
      p_value = test$p.value
    )
  }))
})) %>%
  mutate(
    p_adjusted_bh = p.adjust(p_value, method = "BH"),
    significant_bh_0_05 = p_adjusted_bh < 0.05
  )

write.csv(
  summary_table,
  "output/score_validation/round4_demographic_metric_summary.csv",
  row.names = FALSE
)

write.csv(
  pairwise_results,
  "output/score_validation/round4_demographic_pairwise_effects.csv",
  row.names = FALSE
)

print(summary_table)
print(pairwise_results)
