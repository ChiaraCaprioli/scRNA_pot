
############ explore explanatory variables by PCs ############

ExplVarByPC <- function(sce, var, ndim, mode) {
  
  if (mode == "unscaled") {
    # run PCA on unscaled matrix 
    sce <- scater::runPCA(
      sce, exprs_values = "logcounts",
      ncomponents = ndim, ntop = 2000, scale = F
    )
  }

  # plot
  p1 <- scater::plotExplanatoryVariables(
    sce,
    exprs_values = "logcounts",
    variables = var
    ) +
    theme_bw() +
    scale_color_manual(name = 'variable', values = colors$discrete) +
    theme(panel.grid = element_blank())
  
  x <- getExplanatoryPCs(sce, dimred = "PCA", n_dimred = ndim, variables = var) %>% as.data.frame()
  x <- x %>%
    rownames_to_column(var = 'PC') 
  x <- x %>% pivot_longer(cols = 2:ncol(x), names_to = 'variable', values_to = 'perc_var_explained') 
  x$PC <- factor(x$PC, levels = unique(x$PC))
  x$variable <- factor(x$variable, levels = var)
  
  p2 <- ggplot(x, aes(x = PC, y = perc_var_explained, group = variable, color = variable)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    xlab(paste0("PC1-", ndim)) +
    scale_x_discrete(breaks = c('PC1', 'PC10', 'PC20', 'PC30', 'PC40', 'PC50')) +
    scale_y_log10() +
    ylab('% variance explained (log10)') +
    scale_color_manual(values = colors$discrete) +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
      )

  L <- list(p1,p2)
  return(L)
  
}

