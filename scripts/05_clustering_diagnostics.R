
################ Clustering diagnostics ################

################ PART1: cluster separation ################

## Description
# For any given k/res combination, we evaluate cluster separation by different methods 
# (i.e., silhouette width, purity and modularity).

# Silhouette width evaluates cluster separation.
# For each cell, we compute the average distance to all cells in the same cluster and the average distance 
# to all cells in another cluster, taking the minimum of the averages across all other clusters. 
# The silhouette width for each cell is then defined as the difference between these two values divided by their maximum.
# As a result, cells with large positive silhouette widths are closer to other cells in the same cluster 
# than to cells in different clusters (i.e., clusters with large positive silhouette widths are well-separated from other clusters).
# Moreover, low widths may still occur in well-separated clusters if the internal heterogeneity is high, 
# possibly indicating underclustering.
# As dealing with large datasets makes the calculation of pairwise distances extremely time-consuming,
# we use an approximate approach (implemented in `bluster`) which computes the root of the average squared distances.

# Purity quantifies the degree to which cells from multiple clusters intermingle in expression space. 
# Differently from silhouette width, purity ignores the intra-cluster variance. 
# For each cell, purity is defined as the proportion of neighboring cells that are assigned to the same cluster, 
# after some weighting to adjust for differences in the number of cells between clusters. 
# Well-separated clusters should exhibit little intermingling and thus high purity values for all member cells.

# WCSS quantifies how close data points within a cluster are to the centroid of that cluster.
# It is computed as the sum of squares of the distances of each data point in all clusters to their respective centroids.
# The lower the sum, the more coherent the cluster (i.e., there is no internal structure to score with subclustering).

# Modularity is defined as the (scaled) difference between the observed total weight of edges between nodes in 
# the same cluster and the expected total weight if edge weights were randomly distributed across all pairs of nodes. 
# Larger modularity values indicate that there most edges occur within clusters, suggesting that the clusters are 
# sufficiently well separated to avoid edges forming between neighboring cells in different clusters.
# We use the pairwiseModularity() function from bluster to get the ratio of the observed to expected 
# sum of weights between each pair of clusters (ratio is less dependent on the number of cells in each cluster).
# The function returns a matrix where each row/column corresponds to a cluster and each entry contains 
# the ratio of the observed to total weight of edges between cells in the respective clusters. 
# Concentration of the weight on the diagonal indicates that most  of the clusters are well-separated, 
# while some modest off-diagonal entries represent closely related clusters with more inter-connecting edges.

#' @param seurat Input seurat object.
#' @param k_list List of k to evaluate.
#' @param res_list List of resolutions to evaluate.
#' @param seed Random number generator.
#' @param path Path where results will be saved.
#' @param ndims Number of dimensions to consider.
#' @param method Method to evaluate cluster separation. Can be one of silhouette, purity, WCSS, modularity or all.
#' @return Plots specific to each method to evaluate each k/res combination.

ClusterDiagnostics1 <- function(seurat, k_list, res_list, seed, path, ndims, metric) {
  
  ###### start logger ######
  start_time <- Sys.time()
  info(my_logger, paste0("Evaluating part 1"))
  
  ###### set seed ######
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ###### silhouette ######
  if (metric == "all" | metric == "silhouette") {
    
    info(my_logger, paste0("Evaluating silhouette..."))
    
    L_mean_sil <- list()
    for (k in k_list) {
      
      L1 <- list()
      for (r in res_list) {
        
        res <- paste0("RNA_snn_k_", k, "_res_", r)
        sil.approx <- bluster::approxSilhouette(
          Embeddings(seurat[['scanorama']]), 
          clusters = seurat@meta.data[[res]]
        )
        sil.data <- as.data.frame(sil.approx)
        sil.data$closest <- factor(
          ifelse(sil.data$width > 0, seurat@meta.data[[res]], sil.data$other)
        )
        sil.data$cluster <- as.factor(seurat@meta.data[[res]])
        
        p1 <- ggplot(sil.data, aes(x = cluster, y = width, group = cluster)) +
          geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
          geom_hline(yintercept = mean(sil.data$width), color = "red", linetype = 2) +
          theme_bw() +
          xlab("cluster_id") +
          scale_y_continuous(limits = c(-1,1)) +
          ggtitle(paste0("k_", k, "_res_", r)) +
          theme(
            panel.grid = element_blank(),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 8)
          )
        
        L1[[paste0("res_", r)]] <- p1
        
        L_mean_sil[[paste0("k_", k, "_res_", r)]] <- mean(sil.data$width)
        
      }
      
      k_plots1 <- plot_grid(plotlist = L1, align = "hv")
      
      ggsave(
        paste0(path, "/plots/k", k, "_silhouette.png"),
        k_plots1 + 
          plot_annotation(
            title = paste0("k = ", k),
            theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5, vjust = 2))
          ),
        width = 12, height = 4*(length(res_list)/2)
      )
      
    }
    
    info(my_logger, paste0("Done!"))
  
  }
  
  ###### purity ######
  if (metric == "all" | metric == "purity") {
    
    info(my_logger, paste0("Evaluating purity..."))
    
    L_mean_purity <- list()
    for (k in k_list) {
      
      L2 <- list()
      for (r in res_list) {
        
        res = paste0("RNA_snn_k_", k, "_res_", r)
        pur_data <- neighborPurity(
          Embeddings(seurat[['scanorama']]), 
          clusters = seurat@meta.data[[res]]
        )
        
        pur_data <- as.data.frame(pur_data)
        pur_data$maximum <- factor(pur_data$maximum)
        pur_data$cluster <- as.factor(seurat@meta.data[[res]])
        
        p2 <- ggplot(pur_data, aes(x = cluster, y = purity, group = cluster)) +
          geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
          geom_hline(yintercept = mean(pur_data$purity), color = "red", linetype = 2) +
          theme_bw() +
          xlab("cluster_id") +
          ggtitle(paste0("k_", k, "_res_", r)) +
          theme(
            panel.grid = element_blank(),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 8)
          )
        
        L2[[paste0("res_", r)]] <- p2
        
        L_mean_purity[[paste0("k_", k, "_res_", r)]] <- mean(pur_data$purity)
        
      }
      
      k_plots2 <- plot_grid(plotlist = L2, align = "hv")
      
      ggsave(
        paste0(path, "/plots/k", k, "_purity.png"),
        k_plots2 + 
          plot_annotation(
            title = paste0("k = ", k),
            theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5, vjust = 2))
          ),
        width = 12, height = 4*(length(res_list)/2)
      )
    }
    info(my_logger, paste0("Done!"))
    
  }
  
  ###### WCSS ######
  if (metric == "all" | metric == "WCSS") {
    
    info(my_logger, paste0("Evaluating within-cluster sum of squares..."))
    
    L_sum_wcss <- list()
    for (k in k_list) {
      
      L3 <- list()
      for (r in res_list) {
        
        res = paste0("RNA_snn_k_", k, "_res_", r)
        
        wcss <- as.data.frame(
          clusterRMSD(
            x = Embeddings(seurat[['scanorama']]), 
            clusters = seurat@meta.data[[res]],
            sum = TRUE
            )
          )
        colnames(wcss) <- "wcss"
        wcss$cluster <- as.character(1:nrow(wcss))
        wcss$cluster <- factor(wcss$cluster, levels = 1:nrow(wcss))
        
        p3 <- ggplot(wcss, aes(x = cluster, y = wcss)) +
          geom_point() +
          geom_line(group = "cluster") +
          geom_hline(yintercept = mean(wcss$wcss), color = "red", linetype = 2) +
          theme_bw() +
          xlab("cluster_id") +
          ylab("WCSS") +
          ggtitle(paste0("k_", k, "_res_", r)) +
          theme(
            panel.grid = element_blank(),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 8)
          )
        
        L3[[paste0("res_", r)]] <- p3
        
        L_sum_wcss[[paste0("k_", k, "_res_", r)]] <- sum(wcss$wcss)
        
      }
      
      k_plots3 <- plot_grid(plotlist = L3, align = "hv")
      
      ggsave(
        paste0(path, "/plots/k", k, "_WCSS.png"),
        k_plots3 + 
          plot_annotation(
            title = paste0("k = ", k),
            theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5, vjust = 2))
          ),
        width = 12, height = 4*(length(res_list)/2)
      )
      
    }
    info(my_logger, paste0("Done!"))
  }
  
  ###### modularity ######
  if (metric == "all" | metric == "modularity") {
    
    info(my_logger, paste0("Evaluating modularity..."))
    
    sce <- as.SingleCellExperiment(seurat)
    reducedDim(sce, 'scanorama') <- Embeddings(seurat[['scanorama']])[, 1:ndims]
    g <- scran::buildSNNGraph(sce, use.dimred = 'scanorama')
    
    for (k in k_list) {
      
      L4 <- list()
      for (r in res_list) {
        res = paste0("RNA_snn_k_", k, "_res_", r)
        ratio <- bluster::pairwiseModularity(g, seurat@meta.data[[res]], as.ratio = TRUE)
        ratio_to_plot <- log10(ratio+1)
        x = ratio_to_plot %>%
          as_tibble() %>%
          rownames_to_column(var = 'cluster_1') %>%
          pivot_longer(
            cols = 2:ncol(.),
            names_to = 'cluster_2',
            values_to = 'probability') %>%
          mutate(
            cluster_1 = factor(cluster_1, levels = rev(unique(cluster_1))),
            cluster_2 = factor(cluster_2, levels = unique(cluster_2))) 
        
        p4 <- ggplot(x, aes(cluster_2, cluster_1, fill = probability)) +
          geom_tile(color = 'white') +
          #geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
          scale_x_discrete(name = 'cluster_id', position = 'top') +
          scale_y_discrete(name = 'cluster_id') +
          ggtitle(paste0("k_", k, "_res_", r)) +
          scale_fill_gradient(
            name = 'log10(ratio)', low = 'white', high = '#c0392b', na.value = '#bdc3c7',
            guide = guide_colorbar(
              frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
              title.theme = element_text(hjust = 1, angle = 90),
              barwidth = 0.75, barheight = 10)
          ) +
          coord_fixed() +
          theme_bw() +
          theme(
            legend.position = 'right',
            panel.grid.major = element_blank(),
            axis.text = element_text(size = 5, color = "black"),
            axis.ticks = element_line(linewidth = 0.5)
          )
        
        L4[[paste0("res_", r)]] <- p4
      }
      
      k_plots4 <- plot_grid(plotlist = L4, align = "hv")
      
      ggsave(
        paste0(path, "/plots/k", k, "_modularity.png"),
        k_plots4 + 
          plot_annotation(
            title = paste0("k = ", k),
            theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5, vjust = 2))
          ),
        width = 12, height = 4*(length(res_list)/2)
      )
    }
    info(my_logger, paste0("Done!"))

  }
  
  ###### all metrics vs n clusters ######
  if (metric == "all") {
    
    info(my_logger, paste0("Evaluating all metrics vs n clusters..."))
    
    cluster_all <- seurat@meta.data %>% 
      dplyr::select(
        contains(paste0("k_",k_list, "_")) &
          contains(paste0("res_",res_list)) 
      )
    n_clusters <- apply(cluster_all,2,max) %>% 
      as.data.frame() %>%
      rownames_to_column(var = "k_res") %>%
      dplyr::rename("n_clusters" = 2)
    n_clusters$k_res <- str_remove(n_clusters$k_res, "RNA_snn_")
    n_clusters$k_res <- factor(n_clusters$k_res, levels = unique(n_clusters$k_res))
    
    mean_sil <- do.call(rbind, L_mean_sil) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "k_res")
    colnames(mean_sil)[2] <- "mean_silhouette"
    
    mean_pur <- do.call(rbind, L_mean_purity) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "k_res")
    colnames(mean_pur)[2] <- "mean_purity"
    
    sum_wcss <- do.call(rbind, L_sum_wcss) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "k_res")
    colnames(sum_wcss)[2] <- "sum_wcss"
    
    df_k <- Reduce(full_join, list(n_clusters, mean_sil, mean_pur, sum_wcss)) %>%
      pivot_longer(cols = 2:ncol(.), names_to = "metric", values_to = "score")
    
    df_k$metric <- str_replace_all(df_k$metric, pattern = "_", replacement = " ")
    df_k$metric <- factor(df_k$metric, levels = unique(df_k$metric))
    df_k$k_res <- factor(df_k$k_res, levels = unique(df_k$k_res))  
    
    write.csv(df_k, paste0(path, "/tables/cluster_separation_metrics.csv"))
    
    p <- df_k %>%
      ggplot(aes(k_res, score, color = metric)) +
      geom_point() +
      geom_line(group = "resolution") +
      ylab('Score') +
      theme_bw() +
      facet_wrap(~metric, scales = "free") +
      theme(
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
        panel.grid = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none"
      )
    
    ggsave(
      paste0(path, "/plots/all_metrics.png"),
      p,
      width = 9, height = 3*(length(res_list)/2)
    )
    
    info(my_logger, paste0("Done!"))
    
  }
  
  ###### end logger ######
  end_time <- Sys.time()
  
  info(my_logger, 
       paste0("Part 1 successfully done in ", 
              floor(as.numeric(end_time - start_time, units = "mins")), " mins")
  )
  
}

################ PART2: evaluate cluster trees across different res within a given k   ################

## Description
# For any given k, we evaluate the relationships between different clustering resolutions. 

# We use the clustree package to capture the redistribution of cells from one clustering to another 
# at progressively higher resolution.
# The size of the nodes is proportional to the number of cells in each cluster, and the edges 
# depict cells in one cluster that are reassigned to another cluster at a different resolution. 
# The color of the edges is defined according to the number of reassigned cells and the opacity is 
# defined from the corresponding proportion relative to the size of the lower-resolution cluster.

ClusterDiagnostics2 <- function(seurat, k_list, seed, path) {
  
  start_time <- Sys.time()
  info(my_logger, paste0("Evaluating part 2"))
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ###### cluster tree ######
  info(my_logger, paste0("Evaluating cluster tree..."))
  
  for (k in k_list) {
    p <- clustree(
      seurat, 
      prefix = paste0("RNA_snn_k_", k, "_res_"),
      node_text_size = 0
      ) +
      scale_color_manual(values = colors$discrete) +
      theme(
        legend.position = "right", 
        legend.text = element_text(size = 20)
        )
    ggsave(paste0(path, "/plots/k", k, "_clustree.png"), p, width = 20, height = 20)
    
  }
  
  info(my_logger, paste0("Done!"))
  
  end_time <- Sys.time()
  info(my_logger, 
       paste0("Part 2 successfully done in ", 
              floor(as.numeric(end_time - start_time, units = "mins")), " mins")
  )
  
}

################ PART3: evaluate all k/res combinations ################
ClusterDiagnostics3 <- function(seurat, k_list, res_list, seed, path, ndims) {
  
  start_time <- Sys.time()
  info(my_logger, paste0("Evaluating part 3"))
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ###### ARI ######
  info(my_logger, paste0("Evaluating ARI...")) 
  
  ### ARI
  cluster_all <- seurat@meta.data %>% 
    dplyr::select(
      contains(paste0("k_",k_list)) &
        contains(paste0("res_",res_list)) 
    )
  
  res <- matrix(nrow = ncol(cluster_all), ncol = ncol(cluster_all))
  
  for (i in (1:ncol(cluster_all))) {
    
    for (j in (i:ncol(cluster_all))){
      
      res[i,j] <- res[j,i] <- pairwiseRand(cluster_all[,i], cluster_all[,j], mode="index", adjusted = TRUE)
      
    }
  }
  
  colnames(res) <- str_remove(colnames(cluster_all), "RNA_snn_") 
  rownames(res) <- str_remove(colnames(cluster_all), "RNA_snn_") 
  
  pdf(paste0(path, "/plots/ARI_heat.pdf"), width = 8, height = 8)
  h = Heatmap(
    res, 
    name = "ARI",
    cluster_rows = T,
    cluster_columns = T,
    show_row_dend = F,
    show_column_dend = F,
    show_row_names = T,
    row_names_side = "left",
    show_column_names = T,
    column_names_side = "top",
    col = viridisLite::mako(n = 50, direction = 1),
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    width = unit(12, "cm"), height = unit(12, "cm")
  )
  print(h)
  dev.off()
  
  info(my_logger, paste0("Done!"))

  end_time <- Sys.time()
  info(my_logger, 
       paste0("Part 3 successfully done in ", 
              floor(as.numeric(end_time - start_time, units = "mins")), " mins")
  )
  
}

################ PART4: relationship between clusters and cell type annotation ################

## Description
# For any given k/res combination, 

ClusterDiagnostics4 <- function(seurat, k_list, res_list, seed, path, var) {
  
  ###### Start logger ######
  start_time <- Sys.time()
  info(my_logger, paste0("Evaluating part 4"))
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ###### Cell type fraction across clusters ######
  info(my_logger, paste0("Evaluating cell type fraction across clusters...")) 
  
  for (k in k_list) {
    
    ht_list = NULL
    for (r in res_list) {
      
      res <- paste0("RNA_snn_k_", k, "_res_", r)
    
      mat_cl_lin <- seurat@meta.data %>% dplyr::select(c(res, var))
      mat_cl_lin <- table(as.character(mat_cl_lin[[res]]), as.character(mat_cl_lin[[var]]))
      mat_cl_lin <- t(mat_cl_lin)/apply(mat_cl_lin,2,sum) # fraction of cells from same lineage across clusters
      
      c_order <- sort(as.numeric(colnames(mat_cl_lin)))
      r_order <- levels(seurat@meta.data[[var]])
      
      mat_cl_lin <- mat_cl_lin[match(r_order, rownames(mat_cl_lin)), match(c_order, colnames(mat_cl_lin))]

      ht_list = ht_list + 
        Heatmap(
          mat_cl_lin, 
          name = "Fraction cell type",
          column_title = str_remove_all(res, pattern = "RNA_snn_"),
          cluster_rows = F,
          cluster_columns = F,
          show_row_dend = F,
          show_column_dend = F,
          show_row_names = T,
          row_names_side = "left",
          show_column_names = T,
          column_names_side = "bottom",
          col = viridisLite::mako(n = 50, direction = -1),
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10),
          column_names_rot = 90,
          width = unit(10, "cm"), height = unit(10, "cm")
          )

  }
  
  pdf(paste0(path, "/plots/k", k, "_cluster_", var, ".pdf"), width = 6*length(res_list), height = 7)
  draw(
    ht_list, 
    column_title = paste0("k = ", k), 
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    ht_gap = unit(0.3, "cm"),
    auto_adjust = FALSE
  )
  dev.off()
  
  }
  
  info(my_logger, paste0("Done!"))
  
  ###### ARI ######
  L_ari <- list()
  for (k in k_list) {
    for (r in res_list) {
      res = paste0("RNA_snn_k_", k, "_res_", r)
      
      df_ARI <- data.frame(
        cluster = unlist(seurat[[res]]),
        lineage = seurat@meta.data[[var]]
      )
      ari <- ARI(df_ARI$cluster, df_ARI$lineage)
      L_ari[[res]] <- ari
    }
  }
  
  ari_df <- do.call(rbind, L_ari) %>%
    as.data.frame() %>%
    rownames_to_column(var = "cluster") %>%
    dplyr::rename("ARI" = "V1")
  
  ari_df$cluster <- str_remove(ari_df$cluster, "RNA_snn_")
  
  p <- ggplot(ari_df, aes(x = reorder(cluster, ARI), y = ARI)) +
    geom_point() +
    geom_line(group = "ARI") +
    theme_bw() +
    ggtitle("Agreement between clustering and cell lineage") +
    theme(panel.grid = element_blank(),
          aspect.ratio = 0.7,
          axis.title.x = element_blank(),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
          plot.title = element_text(face = "bold", hjust = 0.5, vjust = 2))
  
  ggsave(paste0(path, "/plots/ari_cluster_", var, ".png"),
         p,
         width = 6, height = 4)
  
  ###### End logger ######
  end_time <- Sys.time()
  info(my_logger, 
       paste0("Part 4 successfully done in ", 
              floor(as.numeric(end_time - start_time, units = "mins")), " mins")
  )
  
}

################ TO DO ################
# stability by bootstrapping
# qc metrics by cluster on chosen solutions
# UMAP on chosen solutions

