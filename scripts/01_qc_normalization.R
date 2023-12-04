
########## Score doublets ##########

run_scDblFinder <- function(seurat) {
  
  # Convert seurat to sce
  sce <- as.SingleCellExperiment(seurat)
  
  # Score doublets independently for each sample
  sce <- scDblFinder(
    sce,
    samples = sce$sample_id,
    multiSampleMode = 'split',
    BPPARAM = BiocParallel::SerialParam()
  )
  
  # Back to seurat
  seurat <- as.Seurat(sce)
  return(seurat)
  
}

########## Plot before/after removing doublets ##########

PlotBeforeAfterDoublets <- function(data) {
  
  # nCount
  p1 <- ggplot(data, aes(x = scDblFinder.class, y = nCount_RNA, fill = scDblFinder.class)) +
    geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
    theme_bw() +
    scale_fill_manual(values = colors$discrete) +
    scale_x_discrete(limits = rev(levels(data[["scDblFinder.class"]]))) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = "Linear scale", y = "N transcripts") +
    theme(
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) 
  
  p2 <- ggplot(data, aes(x = scDblFinder.class, y = nCount_RNA, fill = scDblFinder.class)) +
    geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
    theme_bw() +
    scale_fill_manual(values = colors$discrete) +
    scale_x_discrete(limits = rev(levels(data[["scDblFinder.class"]]))) +
    scale_y_log10(labels = scales::comma) +
    labs(title = "Log scale", y = "N transcripts") +
    theme(
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) 
  
  p_count <- p1 + p2 +
  plot_annotation(
    title = 'N transcripts',
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
  )
  
  # nFeature
  p3 <- ggplot(data, aes(x = scDblFinder.class, y = nFeature_RNA, fill = scDblFinder.class)) +
    geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
    theme_bw() +
    scale_fill_manual(values = colors$discrete) +
    scale_x_discrete(limits = rev(levels(data[["scDblFinder.class"]]))) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = "Linear scale", y = "N genes") +
    theme(
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) 
  
  p4 <- ggplot(data, aes(x = scDblFinder.class, y = nFeature_RNA, fill = scDblFinder.class)) +
    geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
    theme_bw() +
    scale_fill_manual(values = colors$discrete) +
    scale_x_discrete(limits = rev(levels(data[["scDblFinder.class"]]))) +
    scale_y_log10(labels = scales::comma) +
    labs(title = "Linear scale", y = "N genes") +
    theme(
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) 
  
  p_feat <- p3 + p4 +
    plot_annotation(
      title = 'N expressed genes',
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
    )
  
  L_plots <- list(p_count, p_feat)
  
  plot_grid <- plot_grid(plotlist = L_plots, ncol = 1)
  
  return(plot_grid)
  
}

########## Adaptive thresholds ##########

# Define thresholds
DefineAdaptiveThresholds <- function(data, min_umi, min_gene, n_mad) {
  
  data <- data %>%
    filter(scDblFinder.class == "singlet")
  
  L_thr <- list()
  for (s in levels(data$sample_id)) {
    
    data_s <- data[which(data$sample_id == s),]
    
    thresholds <- as.data.frame(
      list(
        sample_id = s,
        thresholds_nCount_min = min_umi, 
        thresholds_nCount_max = median(data_s$nCount_RNA) + n_mad*mad(data_s$nCount_RNA),
        thresholds_nFeature_min = min_gene,
        thresholds_nFeature_max = median(data_s$nFeature_RNA) + n_mad*mad(data_s$nFeature_RNA),
        thresholds_percent_MT_min = 0, 
        thresholds_percent_MT_max = 25
      )
    ) 
    
    L_thr[[s]] <- thresholds
    
  }
  
  thresholds <- do.call(rbind, L_thr)
  return(thresholds)
  
}

# Diagnostic plots
PlotAdaptiveThresholds <- function(data, thresholds) {
  
  data <- data %>%
    filter(scDblFinder.class == "singlet")
  
  ## UMI over genes and % MT 
  L_plots <- list()
  for (s in levels(data$sample_id)) {
    
    thresholds_s <- thresholds[which(thresholds$sample_id == s),]
    
    p_s <- data %>%
      filter(sample_id == s &
               percent_MT < thresholds_s$thresholds_percent_MT_max) %>%
      ggplot(aes(nCount_RNA, nFeature_RNA, color = percent_MT)) +
      geom_point(size = 0.1) +
      scale_x_continuous(name = "N transcripts", labels = scales::comma) +
      scale_y_continuous(name = "N expressed genes", labels = scales::comma) +
      geom_vline(xintercept = thresholds_s$thresholds_nCount_min, color = "red") +
      geom_vline(xintercept = thresholds_s$thresholds_nCount_max, color = "red") +
      geom_hline(yintercept = thresholds_s$thresholds_nFeature_min, color = "red") +
      geom_hline(yintercept = thresholds_s$thresholds_nFeature_max, color = "red") +
      theme_bw() +
      ggtitle(s) +
      scale_color_viridis(
        name = "% MT transcripts",
        limits = c(0,100),
        labels = scales::label_percent(scale = 1),
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
      ) +
      theme(
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", hjust = 0.5, vjust = 2),
        aspect.ratio = 1
      )
    
    L_plots[[s]] <- p_s
    
  }
  
  plot_grid <- plot_grid(plotlist = L_plots, ncol = 3)
  ggsave(paste0(PATH_RESULTS, "/plots/adaptive_thresholds_umi_over_feat_mt.png"), plot_grid, width = 15, height = 15)
  
  ## Violin plots
  
  # prepare data
  data <- full_join(data, thresholds)
  data$sample_id <- factor(
    data$sample_id, 
    levels = data[order(data$cohort), ] %>%
      pull(sample_id) %>%
      unique()
  )
  
  # nCounts
  p1 <- ggplot(data, aes(x = sample_id, y = nCount_RNA, fill = sample_id)) + 
    facet_wrap(~factor(sample_id), scales = "free_x", nrow = 1) +
    geom_violin(draw_quantiles = c(0.5), scale = "area", trim = FALSE) +
    theme_bw() +
    scale_fill_manual(values = colors$sample) +
    scale_y_continuous(name = "N transcripts", labels = scales::comma) +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      strip.text = element_text(face = "bold")
    ) +
    geom_hline(
      data = thresholds,
      aes(yintercept = thresholds_nCount_min, color="red")
    ) +
    geom_hline(
      data = thresholds, 
      aes(yintercept = thresholds_nCount_max, color="red")
    ) 
  
  p2 <- ggplot(data, aes(x = sample_id, y = nCount_RNA, fill = sample_id)) + 
    facet_wrap(~factor(sample_id), scales = "free_x", nrow = 1) +
    geom_violin(draw_quantiles = c(0.5), scale = "area", trim = FALSE) +
    theme_bw() +
    scale_fill_manual(values = colors$sample) +
    scale_y_log10(name = "N transcripts (log10)", labels = scales::comma) +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      strip.text = element_text(face = "bold")
    ) + 
    geom_hline(
      data = thresholds,
      aes(yintercept = thresholds_nCount_min, color="red")
    ) +
    geom_hline(
      data = thresholds, 
      aes(yintercept = thresholds_nCount_max, color="red")
    ) 
  
  vln_count <- p1 / p2 +
    plot_annotation(
      title = "N transcripts",
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, vjust = 2)
        )
    )
  ggsave(paste0(PATH_RESULTS, "/plots/vln_adaptive_thresholds_UMI.png"), 
         vln_count, width = length(unique(data$sample_id)), height = 4)
  
  # nFeature
  p3 <- ggplot(data, aes(x = sample_id, y = nFeature_RNA, fill = sample_id)) + 
    facet_wrap(~factor(sample_id), scales = "free_x", nrow = 1) +
    geom_violin(draw_quantiles = c(0.5), scale = "area", trim = FALSE) +
    theme_bw() +
    scale_fill_manual(values = colors$sample) +
    scale_y_continuous(name = "N genes", labels = scales::comma) +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      strip.text = element_text(face = "bold")
    ) +
    geom_hline(
      data = thresholds,
      aes(yintercept = thresholds_nFeature_min, color="red")
    ) +
    geom_hline(
      data = thresholds, 
      aes(yintercept = thresholds_nFeature_max, color="red")
    ) 
  
  p4 <- ggplot(data, aes(x = sample_id, y = nFeature_RNA, fill = sample_id)) + 
    facet_wrap(~factor(sample_id), scales = "free_x", nrow = 1) +
    geom_violin(draw_quantiles = c(0.5), scale = "area", trim = FALSE) +
    theme_bw() +
    scale_fill_manual(values = colors$sample) +
    scale_y_log10(name = "N genes (log10)", labels = scales::comma) +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      strip.text = element_text(face = "bold")
    ) + 
    geom_hline(
      data = thresholds,
      aes(yintercept = thresholds_nFeature_min, color="red")
    ) +
    geom_hline(
      data = thresholds, 
      aes(yintercept = thresholds_nFeature_max, color="red")
    ) 
  
  vln_feat <- p3 / p4 +
    plot_annotation(
      title = "N genes",
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, vjust = 2)
    )
    )
  ggsave(paste0(PATH_RESULTS, "/plots/vln_adaptive_thresholds_features.png"), 
         vln_feat, width = length(unique(data$sample_id)), height = 4)
  
  # % MT genes
  p5 <- ggplot(data, aes(x = sample_id, y = percent_MT, fill = sample_id)) + 
    facet_wrap(~factor(sample_id), scales = "free_x", nrow = 1) +
    geom_violin(draw_quantiles = c(0.5), scale = "area", trim = FALSE) +
    theme_bw() +
    scale_fill_manual(values = colors$sample) +
    scale_y_continuous(
      name = "% MT transcripts",
      limits = c(0,100),
      labels = scales::label_percent(scale = 1)
    ) +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      strip.text = element_text(face = "bold")
    ) +
    geom_hline(
      data = thresholds,
      aes(yintercept = thresholds_percent_MT_min, color="red")
    ) +
    geom_hline(
      data = thresholds, 
      aes(yintercept = thresholds_percent_MT_max, color="red")
    ) 
  
  vln_mt <- p5 +
    plot_annotation(
      title = "% MT transcripts",
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, vjust = 2)
        )
    )
  ggsave(paste0(PATH_RESULTS, "/plots/vln_adaptive_thresholds_mt.png"), 
         vln_mt, width = length(unique(data$sample_id)), height = 2)
  
}

# Check QC dropout
CheckQCDropout <- function(data, thresholds) {
  
  L_filtered_cells <- list()
  L_cells_qc <- list()
  for (s in levels(data$sample_id)) {
    
    thresholds_s <- thresholds[which(thresholds$sample_id == s),]
    
    data_s <- data[which(data$sample_id == s),]
    data_s_after <- data_s %>%  
      filter(
        scDblFinder.class == "singlet",
        nCount_RNA >= thresholds_s$thresholds_nCount_min,
        nCount_RNA <= thresholds_s$thresholds_nCount_max,
        nFeature_RNA >= thresholds_s$thresholds_nFeature_min,
        nFeature_RNA <= thresholds_s$thresholds_nFeature_max,
        percent_MT >= thresholds_s$thresholds_percent_MT_min,
        percent_MT <= thresholds_s$thresholds_percent_MT_max
      )
    L_filtered_cells[[s]] <- data_s_after    
    
    cells_qc <- tibble(
      sample_id = s,
      before_qc = nrow(data_s),
      after_doublets_removal = nrow(data_s[which(data_s$scDblFinder.class == "singlet"),]),
      after_qc = nrow(data_s_after)
    )
    L_cells_qc[[s]] <- cells_qc
    
  }

  table_cells_qc <- do.call(rbind, L_cells_qc)
  return(table_cells_qc)
  
}

# Apply thresholds
ApplyAdaptiveThresholds <- function(data, thresholds) {
  
  L_flagged_cells <- list()
  L_cells_qc <- list()
  
  for (s in levels(data$sample_id)) {
    
    thresholds_s <- thresholds[which(thresholds$sample_id == s),]
    data_s <- data[which(data$sample_id == s),]
    
    data_s$pass_qc <- ifelse(
      data_s$scDblFinder.class == "singlet" &
      data_s$nCount_RNA >= thresholds_s$thresholds_nCount_min &
      data_s$nCount_RNA <= thresholds_s$thresholds_nCount_max &
      data_s$nFeature_RNA >= thresholds_s$thresholds_nFeature_min &
      data_s$nFeature_RNA <= thresholds_s$thresholds_nFeature_max &
      data_s$percent_MT >= thresholds_s$thresholds_percent_MT_min &
      data_s$percent_MT <= thresholds_s$thresholds_percent_MT_max,
      TRUE, FALSE
    )

    L_flagged_cells[[s]] <- data_s

    cells_qc <- tibble(
      sample_id = s,
      before_qc = nrow(data_s),
      after_doublets_removal = nrow(data_s[which(data_s$scDblFinder.class == "singlet"),]),
      after_qc = nrow(data_s[which(data_s$pass_qc == TRUE),])
    )
    L_cells_qc[[s]] <- cells_qc
    
  }
  
  original_flagged <- do.call(rbind, L_flagged_cells)
  rownames(original_flagged) <- str_split(rownames(original_flagged), "\\.", simplify=T)[,2]
  
  table_cells_qc <- do.call(rbind, L_cells_qc)
  
  L <- list(original_flagged = original_flagged, table_cells_qc = table_cells_qc)
  return(L)
  
}

########## Remove lowly expressed genes ##########

RemoveLowlyExpressedGenes <- function(seurat) {
  
  counts <- seurat@assays$RNA@counts
  genes <- tibble(
    gene = rownames(counts),
    count = Matrix::rowSums(counts),
    cells = Matrix::rowSums(counts != 0)
  )
  genes_to_keep <- genes %>% dplyr::filter(cells >= floor(ncol(seurat)*0.0005)) %>% pull(gene)
  
  message("Filtering ", length(genes_to_keep), " out of ", length(genes$gene), " genes")
  
  seurat <- seurat[genes_to_keep,]
  return(seurat)
  
}


