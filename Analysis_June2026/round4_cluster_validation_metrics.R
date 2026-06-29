local_lib <- file.path(getwd(), "libs", paste0("R-", getRversion()[, 1], ".", getRversion()[, 2]))
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(cluster))

source("0.utils.R")

dir.create("output/score_validation", recursive = TRUE, showWarnings = FALSE)
dir.create("figs", recursive = TRUE, showWarnings = FALSE)

require_analysis_files(c("data_filtered/meta_filtered.csv"))

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
meta_dat <- meta_dat %>% distinct(id, timepoint, .keep_all = TRUE)

cluster_obj <- cluster_om_trajectories(meta_dat, "OMscore_calib", k = 3, grid_step = 0.5)
X2clust <- cluster_obj$dist_matrix
d <- dist(X2clust)
hc <- hclust(d, method = "ward.D2")

k_values <- 2:5

set.seed(8790)
gap <- cluster::clusGap(
  X2clust,
  FUNcluster = function(x, k) {
    hc_tmp <- hclust(dist(x), method = "ward.D2")
    list(cluster = cutree(hc_tmp, k = k))
  },
  K.max = max(k_values),
  B = 100
)
gap_df <- as.data.frame(gap$Tab)
gap_df$k <- seq_len(nrow(gap_df))

choose2 <- function(z) z * (z - 1) / 2

metric_df <- bind_rows(lapply(k_values, function(k) {
  cl <- cutree(hc, k = k)
  sil <- mean(cluster::silhouette(cl, d)[, "sil_width"])
  wss <- sum(vapply(split(seq_len(nrow(X2clust)), cl), function(idx) {
    if (length(idx) <= 1) return(0)
    sum(scale(X2clust[idx, , drop = FALSE], scale = FALSE)^2, na.rm = TRUE)
  }, numeric(1)))
  data.frame(
    k = k,
    gap = gap_df$gap[gap_df$k == k],
    gap_se = gap_df$SE.sim[gap_df$k == k],
    average_silhouette_width = sil,
    total_wss = wss,
    cluster_sizes = paste(as.integer(table(cl)), collapse = "/")
  )
}))

write.csv(
  metric_df,
  "output/score_validation/round4_cluster_validation_metrics.csv",
  row.names = FALSE
)

clusters_3 <- cutree(hc, k = 3)
clusters_4 <- cutree(hc, k = 4)
three_vs_four <- as.data.frame.matrix(table(k3 = clusters_3, k4 = clusters_4))
three_vs_four <- tibble::rownames_to_column(three_vs_four, "k3_cluster")
write.csv(
  three_vs_four,
  "output/score_validation/round4_cluster_validation_k3_vs_k4.csv",
  row.names = FALSE
)

median_final_time <- meta_dat %>%
  group_by(id) %>%
  summarise(final_time = max(timepoint, na.rm = TRUE), .groups = "drop") %>%
  summarise(median_final_time = median(final_time, na.rm = TRUE)) %>%
  pull(median_final_time)

cut_time <- 41
cut_idx <- which(cluster_obj$time_grid == cut_time)
X_cut <- cluster_obj$interpolated[, seq_len(cut_idx), drop = FALSE]
dist_cut <- distMiss(X_cut)
X2clust_cut <- as.matrix(dist_cut)
clusters_cut <- cutree(hclust(dist(X2clust_cut), method = "ward.D2"), k = 3)

cut_compare <- tibble::tibble(
  id = names(clusters_3),
  original_cluster = as.factor(clusters_3),
  restricted_41_day_cluster = as.factor(clusters_cut[names(clusters_3)]),
  same_cluster = original_cluster == restricted_41_day_cluster
)

cut_contingency <- as.data.frame.matrix(
  table(
    original_cluster = cut_compare$original_cluster,
    restricted_41_day_cluster = cut_compare$restricted_41_day_cluster
  )
)
cut_contingency <- tibble::rownames_to_column(cut_contingency, "original_cluster")

write.csv(
  cut_compare,
  "output/score_validation/round4_cluster_validation_41day_assignments.csv",
  row.names = FALSE
)
write.csv(
  cut_contingency,
  "output/score_validation/round4_cluster_validation_41day_contingency.csv",
  row.names = FALSE
)

stability_summary <- data.frame(
  n_patients = nrow(cut_compare),
  median_final_time = median_final_time,
  restricted_time = cut_time,
  n_same_cluster = sum(cut_compare$same_cluster),
  pct_same_cluster = mean(cut_compare$same_cluster) * 100,
  adjusted_rand_index = adjusted_rand_index(
    cut_compare$original_cluster,
    cut_compare$restricted_41_day_cluster
  )
)

write.csv(
  stability_summary,
  "output/score_validation/round4_cluster_validation_41day_summary.csv",
  row.names = FALSE
)

png("figs/round4_cluster_dendrogram.png", width = 1800, height = 900, res = 180)
plot(
  as.dendrogram(hc),
  main = "Trajectory-clustering dendrogram",
  xlab = "Patients",
  ylab = "Height",
  leaflab = "none"
)
rect.hclust(hc, k = 3, border = 2:4)
dev.off()

print(metric_df)
print(three_vs_four)
print(stability_summary)
print(cut_contingency)
