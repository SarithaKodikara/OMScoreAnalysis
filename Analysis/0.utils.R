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