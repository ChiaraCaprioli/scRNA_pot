################ Differential abundance ################

################ Milo ################

# to do: implement design with multiple covariates

RunMilo <- function(
    
  seurat = seurat, 
  contrast = contrast, 
  var = var,
  cov = NULL,
  reduction = reduction, 
  ndims = ndims, 
  alpha = alpha,
  lineage = lineage,
  frequent_lineage = frequent_lineage,
  seed = seed,
  path_save = path_save
  
  ) {
  
  #### Getting prepared ####
  message(paste0(format(Sys.time()), " Getting prepared..."))
  
  # path to save results
  if (!dir.exists(paste0(path_save, "/milo"))) {
    dir.create(paste0(path_save, "/milo"))
  }
  path_milo <- paste0(path_save, "/milo/")
  
  if(!dir.exists(paste0(path_milo, contrast))) {
    dir.create(path = paste0(path_milo, contrast))
  }
  path_save <- paste0(path_milo, contrast)
  
  # set seed
  if (!is.null(seed)) {
    set.seed(seed)
  } else {
    message("Please set random seed")
  }
  
  # subset seurat to contrast and cell types of interest
  obj <- seurat[,which(
    seurat@meta.data[[var]] %in% c(str_split_1(contrast, pattern = "_")[1], str_split_1(contrast, pattern = "_")[2]) &
      seurat@meta.data[[frequent_lineage]] == TRUE
  )]
  obj@meta.data[[var]] <- factor(obj@meta.data[[var]], levels = unique(obj@meta.data[[var]]))
  obj$sample_id <- factor(obj$sample_id, levels = unique(obj$sample_id))
  
  # convert seurat to sce
  sce <- as.SingleCellExperiment(obj, assay = 'RNA')
  
  # create milo object
  data_milo <- Milo(sce)
  
  #### Creating de novo KNN graph ####
  message(paste0(format(Sys.time()), " Creating de novo KNN graph..."))
  
  ## From Milo developers, it is recommended to choose k such that k≥Sx5, 
  ## where S is the number of experimental samples, or to have a distribution 
  ## peaking between 50 and 100, or mean average > 5 x n samples
  
  k = 5*length(unique( data_milo@colData$sample_id ))
  data_milo <- buildGraph(data_milo, k = k, d = ndims, reduced.dim = reduction)
  
  # Define representative neighbourhoods on the KNN graph
  data_milo <- makeNhoods(
    data_milo, 
    prop = 0.1, 
    k = k, d = ndims, 
    refined = TRUE, 
    refinement_scheme = "graph",
    reduced_dims = reduction
    )
  
  #### Checking neighbourhood sizes ####
  message(paste0(format(Sys.time()), " Checking KNN graph..."))
  
  ## This heuristics serves to evaluate whether the value of k used for graph building is appropriate, 
  ## as neighbourhood sizes affect the power of DA testing.
  data_milo <- countCells(
    data_milo, 
    meta.data = as.data.frame(colData(data_milo)), 
    sample = "sample_id"
    )
  
  test <- as.data.frame(data_milo@nhoodCounts) 
  test$size <- rowSums(test)
  
  p = ggplot(test, aes(size)) +
    geom_histogram(bins = 50) +
    theme_bw() +
    xlab("n cells in neighbourhoods") +
    ggtitle(paste0("k = ", k)) +
    geom_vline(xintercept = mean(test$size), linetype = 2, color = 'red') +
    geom_vline(xintercept = 5*length(unique( data_milo@colData$sample_id )), linetype = 2, color = 'orange') +
    scale_x_continuous(breaks = c(0,50,100,max(test$size))) +
    annotate('rect', xmin = 50, xmax = 100, 
             ymin = 0, ymax = Inf, 
             alpha = 0.2, fill = 'lightblue') +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave(paste0(path_save, "/milo_NhoodSizeHist.png"), p, width = 5, height = 5)
  
  #### Setting DA testing ####
  message(paste0(format(Sys.time()), " Setting DA testing..."))

  # prepare design matrix 
  data_design <- data.frame(colData(data_milo))[,c("sample_id", "cohort", "age", "sex")] 
  data_design <- distinct(data_design)
  rownames(data_design) <- data_design$sample_id
  data_design <- data_design[match(colnames(nhoodCounts(data_milo)), rownames(data_design)),]
  data_design[[var]] <- factor(
    data_design[[var]], 
    levels = c(str_split_1(contrast, pattern = "_")[1],str_split_1(contrast, pattern = "_")[2])
  )
  data_design$sample_id <- as.factor(data_design$sample_id)
  
  # set contrast
  contr <- paste0(
    colnames(data_design)[2],str_split_1(contrast, pattern = "_")[1],
    "-",
    colnames(data_design)[2],str_split_1(contrast, pattern = "_")[2]
  )
  
  #### Testing DA ####
  message(paste0(format(Sys.time()), " Testing DA..."))
  
  if (is.null(cov)) {
    
    da_results <- testNhoods(
      data_milo,
      design = ~ 0 + cohort, 
      design.df = data_design, 
      model.contrasts = contr,
      fdr.weighting = "graph-overlap" # https://github.com/MarioniLab/miloR/issues/99
    )
    
  } else {
    
    da_results <- testNhoods(
      data_milo,
      design = ~ 0 + cohort + age, 
      design.df = data_design, 
      model.contrasts = contr,
      fdr.weighting = "graph-overlap"
    )
    
  }
  
  #### Checking DA results ####
  message(paste0(format(Sys.time()), " Checking DA results..."))
  
  if (length(which(da_results$SpatialFDR < alpha)) == 0) {
    
      message(paste0("No nhoods with FDR < ", alpha, " for contrast ", contrast))
    
    } else {

    # Plot neighbourhood graph
    data_milo <- buildNhoodGraph(data_milo)
    
    p <- 
      plotNhoodGraphDA(
        data_milo, 
        da_results, 
        layout = "UMAP",
        alpha = alpha
      ) +
      labs(
        title = paste0(
          str_split_1(contrast, pattern = "_")[1]," Vs ",str_split_1(contrast, pattern = "_")[2]
        ), 
        subtitle = paste0("FDR < ", alpha)) +
      scale_fill_gradient2(name = "logFC",
                           low = "steelblue",
                           mid = "white",
                           high = "darkred",
                           midpoint = 0,
                           space = "Lab",
                           na.value = "lightgrey") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 2),
            plot.subtitle = element_text(face = "italic", hjust = 0.5, vjust = 1))
    
    ggsave(paste0(path_save, "/da_milo_nhood.png"), p,
           width = 7, height = 7)
    
    # check p values
    p1 <- ggplot(da_results, aes(PValue)) + 
      geom_histogram(bins = floor(nrow(da_results)/10)) +
      theme_bw() +
      theme(panel.grid = element_blank())
    
    p2 <- ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
      geom_point(alpha = 0.5, size = 0.5) +
      geom_hline(yintercept = -log10(alpha), color = "red", linetype = 2) +
      theme_bw() +
      theme(panel.grid = element_blank())
    
    ggsave(paste0(path_save, "/da_milo_sig.png"), 
           p1 + p2 +
             patchwork::plot_layout(ncol = 2, widths = 5, heights = 4),
           width = 10, height = 4)
    
    # Find the most abundant cell type within cells in each neighborhood
    da_results <- annotateNhoods(data_milo, da_results, coldata_col = lineage)
    
    p = ggplot(da_results, aes(x = da_results[[9]])) + 
      geom_histogram(bins = floor(nrow(da_results)/10)) +
      geom_vline(xintercept = 0.8, linetype = 2, color = 'red') +
      theme_bw() +
      xlab("lineage") +
      theme(panel.grid = element_blank())
    
    ggsave(paste0(path_save, "/da_milo_lineage_fraction.png"), p,
           width = 5, height = 5)
    
    da_results$lineage <- ifelse(da_results[[9]] < 0.8, "Mixed", da_results[[8]])
    
    plot_data <- da_results %>%
      filter(lineage != "Mixed")
    
    plot_data$lineage <- factor(
      plot_data$lineage,
      levels = levels(obj@meta.data[[lineage]])[which(levels(obj@meta.data[[lineage]]) %in% unique(plot_data$lineage))]
    )
    
    p <- plotDAbeeswarm(plot_data, group.by = "lineage", alpha = alpha) +
      labs(
        title = "Abundance of cell-type neighbourhood",
        subtitle = paste0("FDR < ", alpha)
      ) +
      scale_colour_gradient2(
        low = "steelblue",
        mid = "white",
        high = "darkred",
        midpoint = 0,
        space = "Lab",
        na.value = "lightgrey"
      ) +
      scale_x_discrete(limits = rev(levels(plot_data$lineage))) +
      coord_cartesian(clip = "off") +
      coord_flip(clip = "off") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, vjust = 20, size = 13),
        plot.subtitle = element_text(face = "italic", hjust = 0.5, vjust = 18, size = 12),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        plot.margin = unit(c(3,0.5,0,0), "cm"),
        panel.border = element_rect(linewidth = 0.5),
        aspect.ratio = 1
      ) +
      annotation_custom(grob = grid::linesGrob(
        arrow = arrow(type="open", ends="both", length = unit(3,"mm")), 
        gp = grid::gpar(col = "black", lwd = 1)
      ), 
      xmin = length(unique(plot_data$lineage))+1.3, 
      xmax = length(unique(plot_data$lineage))+1.3, 
      ymin = (min(plot_data$logFC)+1), ymax = (max(plot_data$logFC)-1)) +
      annotation_custom(
        grob = grid::textGrob(
          label = str_split_1(contrast, pattern = "_")[1], 
          hjust = 0.5, vjust = -0.5, 
          rot = 0, gp = grid::gpar(col = "black")
        ), 
        xmin = length(unique(plot_data$lineage))+1.7, 
        xmax = length(unique(plot_data$lineage))+1.7, 
        ymin = 6, ymax = 10) +
      annotation_custom(
        grob = grid::textGrob(
          label = str_split_1(contrast, pattern = "_")[2],
          hjust = 0.5, vjust = -0.5, 
          rot = 0, gp = grid::gpar(col = "black")
        ), 
        xmin = length(unique(plot_data$lineage))+1.7, 
        xmax = length(unique(plot_data$lineage))+1.7, 
        ymin = -10, ymax = -6) 
    
    ggsave(paste0(path_save, "/da_results.png"), p, width = 5, height = 5)
    
  }
  
}


################ Plot scCODA results ################
PlotScCODA <- function(contrast, path_res, alpha, width, height) {
  
  # prepare data
  res <- read_csv(paste0(path_res, contrast, "/", contrast, ".csv"), show_col_types = F)
  res$`Cell Type` <- factor(res$`Cell Type`, levels = res$`Cell Type`)
  res$flag <- case_when(
    res$`Final Parameter` == 0 ~ "ns",
    res$`Final Parameter` != 0 & res$`log2-fold change` > 0 ~ "up",
    res$`Final Parameter` != 0 & res$`log2-fold change` < 0 ~ "down"
  )
  res$flag <- factor(res$flag, levels = c("ns", "up", "down"))
  
  # plot
  p <- ggplot(res, aes(x = res$`Cell Type`, y = res$`log2-fold change`, fill = flag)) +
    geom_col(color = "black", linewidth = 0.3, width = 0.8) +
    theme_bw() +
    scale_fill_manual(values = c("lightgrey", "darkred", "steelblue")) +
    labs(
      title = contrast,
      subtitle = paste0("FDR < ", alpha),
      y = "Log2FC"
    ) +
    coord_cartesian(clip = "off") +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, vjust = 20, size = 13),
      plot.subtitle = element_text(face = "italic", hjust = 0.5, vjust = 18, size = 12),
      axis.title.x = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "none",
      plot.margin = unit(c(2,0,0,3), "cm"),
      panel.border = element_rect(linewidth = 0.5),
      aspect.ratio = 0.5
    ) +
    annotation_custom(
      grob = grid::linesGrob(
        arrow = arrow(type="open", ends="both", length = unit(3,"mm")), 
        gp = grid::gpar(col = "black", lwd = 1)
      ), 
      xmin = -1.8, xmax = -1.8, 
      ymin = min(res$`log2-fold change`), ymax = max(res$`log2-fold change`)
    ) +
    annotation_custom(
      grob = grid::textGrob(
        label = str_split_1(contrast, pattern = "_")[1], 
        hjust = 0.5, vjust = -0.5, just = "top",
        rot = 90, gp = grid::gpar(col = "black", fontsize = 11)
      ), 
      xmin = -2.1, xmax = -2.1, 
      ymin = max(res$`log2-fold change`)-1,
      ymax = max(res$`log2-fold change`)
    ) +
    annotation_custom(
      grob = grid::textGrob(
        label = str_split_1(contrast, pattern = "_")[2],
        hjust = 0.5, vjust = -0.5, just = "bottom",
        rot = 90, gp = grid::gpar(col = "black", fontsize = 11)
      ), 
      xmin = -2.1, xmax = -2.1, 
      ymin = min(res$`log2-fold change`), 
      ymax = min(res$`log2-fold change`)+1
    ) 
  
  # save
  ggsave(paste0(path_res, contrast, "/", contrast, "_da.png"), width = width, height = height)
  
}
