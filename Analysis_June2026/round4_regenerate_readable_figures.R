local_lib <- file.path(getwd(), "libs", paste0("R-", getRversion()[, 1], ".", getRversion()[, 2]))
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages(library(PhiSpace))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(nlme))

source("0.utils.R")

dir.create("figs", recursive = TRUE, showWarnings = FALSE)
panel_dir <- "figure_panels"
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(plot, filename, width, height, dpi = 300) {
  png_path <- file.path(panel_dir, paste0(filename, ".png"))
  pdf_path <- file.path(panel_dir, paste0(filename, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, bg = "white")
}

require_analysis_files(c(
  "data_filtered/meta_filtered.csv",
  "data/taxonomy_table.csv",
  "data_filtered/sce_all.rds",
  "data_filtered/selectedOTUs.rds",
  "output/score_validation/round4_cluster_assignments_from_render.csv",
  "output/score_validation/round4_demographic_merged_readr.csv",
  "output/score_validation/round4_selected_taxa_direction_table.csv",
  "output/score_validation/round4_selected_taxa_supplementary_table.csv"
))

clean_taxonomy <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unclassified"
  x
}

to_numeric <- function(x) {
  as.numeric(gsub("[^0-9.\\-]", "", as.character(x)))
}

direction_label <- function(x) {
  case_when(
    is.na(x) ~ "Not selected",
    x == 1 ~ "Selected positive",
    x == 0 ~ "Selected negative",
    TRUE ~ "Not selected"
  )
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
    coefficient = coef_vec
  )
}

make_loading_plot <- function(coef_df, selected_taxa, title, n_total = 16) {
  selected_taxa <- selected_taxa[!is.na(selected_taxa)]
  taxa_to_show <- coef_df %>%
    arrange(desc(abs(coefficient))) %>%
    pull(taxon)
  taxa_to_show <- unique(c(selected_taxa, taxa_to_show))
  taxa_to_show <- taxa_to_show[seq_len(min(length(taxa_to_show), n_total))]
  coef_df %>%
    filter(taxon %in% taxa_to_show) %>%
    mutate(
      selection = case_when(
        taxon %in% selected_taxa & coefficient >= 0 ~ "Selected positive",
        taxon %in% selected_taxa & coefficient < 0 ~ "Selected negative",
        TRUE ~ "Not selected"
      ),
      taxon = factor(taxon, levels = taxon[order(coefficient)])
    ) %>%
    ggplot(aes(coefficient, taxon, fill = selection)) +
    geom_col(width = 0.78) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey35") +
    scale_fill_manual(
      values = c(
        "Selected positive" = "#D95F5F",
        "Selected negative" = "#4F6BD7",
        "Not selected" = "grey84"
      ),
      limits = c("Selected positive", "Selected negative", "Not selected"),
      breaks = c("Selected positive", "Selected negative", "Not selected")
    ) +
    labs(title = title, x = "PLS coefficient", y = NULL, fill = NULL) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title.x = element_text(size = 13),
      axis.text.y = element_text(size = 10.5),
      axis.text.x = element_text(size = 11),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 11)
    )
}

meta_dat_raw <- read.csv("data_filtered/meta_filtered.csv", row.names = 1)
meta_dat <- meta_dat_raw %>%
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
    heightnew = to_numeric(heightnew),
    weight = to_numeric(weight)
  )
meta_dat <- meta_dat %>% left_join(demographics, by = "id")
rownames(meta_dat) <- rownames(meta_dat_raw)

taxonomy_table <- read.csv("data/taxonomy_table.csv", row.names = 1) %>%
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

direction_table <- read.csv(
  "output/score_validation/round4_selected_taxa_direction_table.csv",
  check.names = FALSE
)
supplementary_table <- read.csv(
  "output/score_validation/round4_selected_taxa_supplementary_table.csv",
  check.names = FALSE
)

global_coef <- get_pls_coefficients(
  microb_mat,
  meta_dat$OMscore_calib,
  taxon_labels
)

cluster_coef <- lapply(1:3, function(cluster_id) {
  cluster_samples <- rownames(meta_dat)[meta_dat$cluster == cluster_id]
  get_pls_coefficients(
    microb_mat[cluster_samples, , drop = FALSE],
    meta_dat[cluster_samples, "OMscore_calib"],
    taxon_labels
  )
})

loading_plots <- list(
  make_loading_plot(
    global_coef,
    direction_table$taxon[!is.na(direction_table$All)],
    "All patients",
    n_total = 18
  ),
  make_loading_plot(
    cluster_coef[[1]],
    direction_table$taxon[!is.na(direction_table$C1)],
    "Cluster 1",
    n_total = 15
  ),
  make_loading_plot(
    cluster_coef[[2]],
    direction_table$taxon[!is.na(direction_table$C2)],
    "Cluster 2",
    n_total = 15
  ),
  make_loading_plot(
    cluster_coef[[3]],
    direction_table$taxon[!is.na(direction_table$C3)],
    "Cluster 3",
    n_total = 15
  )
)

fig3 <- wrap_plots(loading_plots, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")

save_panel(loading_plots[[1]], "fig3a_loading_all_patients", width = 7.2, height = 6.8)
save_panel(loading_plots[[2]], "fig3b_loading_cluster1", width = 7.2, height = 6.4)
save_panel(loading_plots[[3]], "fig3c_loading_cluster2", width = 7.2, height = 6.4)
save_panel(loading_plots[[4]], "fig3d_loading_cluster3", width = 7.2, height = 6.4)

ggsave(
  "figs/round4_Results_fig3.png",
  fig3,
  width = 13.6,
  height = 6.2,
  dpi = 300,
  bg = "white"
)
occurrence_long <- direction_table %>%
  pivot_longer(
    cols = c(All, C1, C2, C3),
    names_to = "context",
    values_to = "direction"
  ) %>%
  mutate(
    context = recode(
      context,
      All = "All",
      C1 = "Cluster 1",
      C2 = "Cluster 2",
      C3 = "Cluster 3"
    ),
    context = factor(context, levels = c("All", "Cluster 1", "Cluster 2", "Cluster 3")),
    association = direction_label(direction),
    taxon = factor(taxon, levels = rev(direction_table$taxon))
  )

occurrence_plot <- ggplot(occurrence_long, aes(context, taxon, fill = association)) +
  geom_tile(color = "grey88", linewidth = 0.35) +
  scale_fill_manual(
    values = c(
      "Selected positive" = "#D95F5F",
      "Selected negative" = "#4F6BD7",
      "Not selected" = "white"
    ),
    limits = c("Selected positive", "Selected negative", "Not selected"),
    breaks = c("Selected positive", "Selected negative", "Not selected")
  ) +
  labs(title = "A. PLSKO-selected taxa", x = NULL, y = NULL, fill = NULL) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 11)
  )

selected_taxa <- direction_table$taxon
selected_otus <- supplementary_table$otu[match(selected_taxa, supplementary_table$taxon)]
names(selected_otus) <- selected_taxa

lmm_meta <- meta_dat %>%
  select(id, timepoint, cluster, age, heightnew, weight) %>%
  mutate(cluster = factor(cluster))
complete_sample_ids <- rownames(lmm_meta)[complete.cases(lmm_meta)]
lmm_meta <- lmm_meta[complete_sample_ids, , drop = FALSE]
lmm_x <- microb_mat[complete_sample_ids, selected_otus, drop = FALSE]
colnames(lmm_x) <- selected_taxa

lmm_terms <- c(
  timepoint = "Timepoint",
  cluster2 = "Cluster 2 vs 1",
  cluster3 = "Cluster 3 vs 1",
  age = "Age",
  heightnew = "Height",
  weight = "Weight"
)

lmm_long <- bind_rows(lapply(selected_taxa, function(taxon) {
  model_df <- data.frame(
    abundance = as.numeric(lmm_x[, taxon]),
    lmm_meta
  )
  fit <- tryCatch(
    nlme::lme(
      abundance ~ timepoint + cluster + age + heightnew + weight,
      random = ~1 | id,
      data = model_df,
      na.action = na.omit,
      control = nlme::lmeControl(opt = "optim")
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(tibble(taxon = taxon, term = names(lmm_terms), beta = NA_real_, p_value = NA_real_))
  }
  tab <- summary(fit)$tTable
  tibble(
    taxon = taxon,
    term = names(lmm_terms),
    beta = tab[names(lmm_terms), "Value"],
    p_value = tab[names(lmm_terms), "p-value"]
  )
})) %>%
  mutate(
    term_label = factor(unname(lmm_terms[term]), levels = unname(lmm_terms)),
    taxon = factor(taxon, levels = rev(selected_taxa)),
    signed_log10_p = sign(beta) * -log10(p_value),
    signed_log10_p = pmax(pmin(signed_log10_p, 12), -12),
    significant = !is.na(p_value) & p_value < 0.05
  )

write.csv(
  lmm_long %>% mutate(taxon = as.character(taxon), term_label = as.character(term_label)),
  "output/score_validation/round4_readable_figure4_lmm_values.csv",
  row.names = FALSE
)

lmm_plot <- ggplot(lmm_long, aes(term_label, taxon)) +
  geom_tile(aes(fill = ifelse(significant, signed_log10_p, NA_real_)),
            color = "grey88", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#B77A20",
    mid = "white",
    high = "#6A3D9A",
    midpoint = 0,
    limits = c(-12, 12),
    na.value = "grey92",
    name = "Significant\nsigned\n-log10(P)"
  ) +
  labs(
    title = "B. LMM predictors of selected-taxon abundance",
    x = NULL,
    y = NULL,
    caption = "Grey cells: P >= 0.05"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.caption = element_text(size = 11, hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 11)
  )

lmm_plot_standalone <- lmm_plot +
  theme(
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_line()
  )

save_panel(occurrence_plot, "fig4a_plsko_selected_taxa_heatmap", width = 6.8, height = 8.2)
save_panel(lmm_plot_standalone, "fig4b_lmm_predictor_heatmap", width = 8.6, height = 8.2)

fig4 <- occurrence_plot + lmm_plot +
  plot_layout(widths = c(1.05, 1.35), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "figs/round4_Results_fig4.png",
  fig4,
  width = 13.2,
  height = 7.6,
  dpi = 300,
  bg = "white"
)
message("Regenerated readable Results_fig3.png and Results_fig4.png plus standalone panels")
