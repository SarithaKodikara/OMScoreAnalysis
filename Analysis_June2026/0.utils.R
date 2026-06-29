# Change gower distance used in cluster::daisy (without standardisation)
distMiss <- function (
    X
) 
{
  X <- data.matrix(X)
  
  Nrows <- nrow(X)
  
  distM <- matrix(NA, nrow = Nrows, ncol = Nrows)
  dimnames(distM) <- list(rownames(X), rownames(X))
  for(irow in 1:(nrow(X)-1)){
    diffMat <- sweep(X[(irow+1):nrow(X),,drop=F], 2, X[irow,], FUN = "-")
    distM[(irow+1):Nrows, irow] <- sqrt(rowMeans(diffMat^2, na.rm = T))
  }
  
  return(as.dist(distM))
  
}

analysis_palette <- function(n) {
  grDevices::hcl.colors(n, palette = "Dark 3")
}

save_result <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("qs", quietly = TRUE)) {
    qs::qsave(x, path)
  } else {
    saveRDS(x, sub("\\.qs$", ".rds", path))
  }
}

read_result <- function(path) {
  if (file.exists(path) && requireNamespace("qs", quietly = TRUE)) {
    qs::qread(path)
  } else {
    rds_path <- sub("\\.qs$", ".rds", path)
    if (!file.exists(rds_path)) {
      stop("No cached result found at ", path, " or ", rds_path, call. = FALSE)
    }
    readRDS(rds_path)
  }
}

result_exists <- function(path) {
  file.exists(path) || file.exists(sub("\\.qs$", ".rds", path))
}

require_analysis_files <- function(paths) {
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0) {
    stop(
      "Missing required analysis input files:\n",
      paste0("  - ", missing_paths, collapse = "\n"),
      "\nRestore the Analysis/data and Analysis/data_filtered folders before rerunning.",
      call. = FALSE
    )
  }
}

plot_nbclust_custom <- function(X, method = c("gap_stat", "silhouette", "wss"),
                                k.max = 10, seed = 8790) {
  method <- match.arg(method)
  k_values <- 1:k.max
  
  if (method == "gap_stat") {
    set.seed(seed)
    gap <- cluster::clusGap(
      X,
      FUNcluster = function(x, k) {
        hc <- hclust(dist(x), method = "ward.D2")
        list(cluster = cutree(hc, k = k))
      },
      K.max = k.max,
      B = 100
    )
    gap_df <- as.data.frame(gap$Tab)
    gap_df$k <- seq_len(nrow(gap_df))
    return(
      ggplot(gap_df, aes(k, gap)) +
        geom_line() +
        geom_point() +
        geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.15) +
        labs(x = "Number of clusters", y = "Gap statistic") +
        theme_bw()
    )
  }
  
  if (method == "silhouette") {
    d <- dist(X)
    vals <- vapply(k_values[-1], function(k) {
      hc <- hclust(d, method = "ward.D2")
      mean(cluster::silhouette(cutree(hc, k = k), d)[, "sil_width"])
    }, numeric(1))
    df <- data.frame(k = k_values[-1], value = vals)
    return(
      ggplot(df, aes(k, value)) +
        geom_line() +
        geom_point() +
        labs(x = "Number of clusters", y = "Average silhouette width") +
        theme_bw()
    )
  }
  
  vals <- vapply(k_values, function(k) {
    hc <- hclust(dist(X), method = "ward.D2")
    cl <- cutree(hc, k = k)
    sum(vapply(split(seq_len(nrow(X)), cl), function(idx) {
      if (length(idx) <= 1) return(0)
      sum(scale(X[idx, , drop = FALSE], scale = FALSE)^2, na.rm = TRUE)
    }, numeric(1)))
  }, numeric(1))
  df <- data.frame(k = k_values, value = vals)
  ggplot(df, aes(k, value)) +
    geom_line() +
    geom_point() +
    labs(x = "Number of clusters", y = "Total within-cluster sum of squares") +
    theme_bw()
}

adjusted_rand_index <- function(x, y) {
  x <- as.factor(x)
  y <- as.factor(y)
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  sum_nij <- sum(choose2(tab))
  sum_ai <- sum(choose2(rowSums(tab)))
  sum_bj <- sum(choose2(colSums(tab)))
  n <- sum(tab)
  total <- choose2(n)
  expected <- sum_ai * sum_bj / total
  max_index <- 0.5 * (sum_ai + sum_bj)
  if (max_index == expected) {
    return(NA_real_)
  }
  (sum_nij - expected) / (max_index - expected)
}

cluster_om_trajectories <- function(meta_dat, score_col, k = 3, grid_step = 0.5) {
  time_grid <- seq(
    min(meta_dat$timepoint, na.rm = TRUE),
    max(meta_dat$timepoint, na.rm = TRUE),
    grid_step
  )
  meta_list <- split(meta_dat, ~id)
  interp_vals <- lapply(
    meta_list,
    function(sub_df) {
      approx(
        x = sub_df$timepoint,
        y = sub_df[[score_col]],
        xout = time_grid
      )
    }
  )
  X_interp <- sapply(interp_vals, function(x) x$y) |> t()
  dist_m <- distMiss(X_interp)
  X2clust <- as.matrix(dist_m)
  hc <- hclust(dist(X2clust), method = "ward.D2")
  clusters <- cutree(hc, k = k)
  list(
    time_grid = time_grid,
    interpolated = X_interp,
    dist = dist_m,
    dist_matrix = X2clust,
    hclust = hc,
    clusters = clusters
  )
}

# Modified loading plot function
loadBarplot_v2 <- function(
    Loadings, comp = "comp1", showInt = F, absVal = T, showNeg = F,
    nfeat = 30, fsize = 14, xlab = "", significant = NULL
){
  
  if(absVal){
    
    plot_dat <- Loadings %>% as.data.frame() %>%
      dplyr::arrange(dplyr::desc(abs(!!sym(comp))))  %>%
      dplyr::slice_head(n = nfeat) %>%
      dplyr::arrange((!!sym(comp)))
  } else {
    
    if(showNeg){
      plot_dat <- Loadings %>% as.data.frame() %>%
        dplyr::arrange(dplyr::desc(-(!!sym(comp))))  %>%
        dplyr::slice_head(n = nfeat) %>%
        dplyr::arrange((!!sym(comp)))
    } else {
      
      plot_dat <- Loadings %>% as.data.frame() %>%
        dplyr::arrange(dplyr::desc(!!sym(comp)))  %>%
        dplyr::slice_head(n = nfeat) %>%
        dplyr::arrange((!!sym(comp)))
    }
  }
  
  if(showInt){
    rname_list <- lapply(
      strsplit(rownames(plot_dat), "<->"),
      function(x) sort(x)
    )
    rname_new <- sapply(
      rname_list,
      function(x) paste(x[1], x[2], sep = "<->")
    )
    rownames(plot_dat) <- rname_new
  }
  
  plot_dat <- plot_dat %>%
    dplyr::mutate(
      interaction = factor(rownames(plot_dat), levels = rownames(plot_dat)),
      fill_color = ifelse((rownames(plot_dat) %in% significant),  !!sym(comp), NA)
    )
  
  if(is.null(significant)){
    p <- plot_dat  %>%
      ggplot() +
      geom_bar(
        aes(x = !!sym(comp), y = interaction, fill = !!sym(comp)),
        stat = "identity"
      ) +
      xlab(xlab) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                           na.value = "gray89") +
      theme_light(base_size = fsize) +
      theme(
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "none"
      )
  }else{
    p <- plot_dat  %>%
      ggplot() +
      geom_bar(
        aes(x = !!sym(comp), y = interaction, fill = fill_color),
        stat = "identity"
      ) +
      xlab(xlab) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                           na.value = "gray89") +
      theme_light(base_size = fsize) +
      theme(
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "none"
      )
  }
  
  return(p)
  
}
