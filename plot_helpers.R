
########## Barplot proportions of grouped variables ##########

# Requires data = table whose columns contain main grouping variable, total cell count and 
# counts by other grouping variables (order is mandatory)

BarGroupProp <- function(data, id.var, var = NULL, col = NULL) {

  data_prop <- data %>%
    dplyr::select(-c('total_cell_count')) %>%
    reshape2::melt(id.vars = id.var) 
  data_prop %>% 
    ggplot(aes(data_prop[[id.var]], value)) +
    geom_bar(
      aes(fill = variable), 
      color = "black", position = 'fill', 
      stat = 'identity', width = 0.7,
      size = 0.3
      ) +
    #geom_text(
    #  data = data,
    #  aes(x = data[[1]], y = Inf, 
    #      label = paste0('n = ', format(total_cell_count, big.mark = ',', trim = TRUE)), 
    #      vjust = 0.5, hjust = -0.2, angle = 90),
    #  color = 'black', size = 3
    #) +
    scale_y_continuous(name = '% cells', labels = scales::percent_format(), expand = c(0.01,0)) +
    scale_fill_manual(name = var, values = col) +
    coord_cartesian(clip = 'off') +
    guides(fill = guide_legend(keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm"))) +
    theme_bw() +
    theme(
      text = element_text(size = 13, color = 'black'),
      aspect.ratio = 1,
      axis.text = element_text(color = 'black'),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      legend.text = element_text(size = 9),
      plot.margin = margin(t = 1, r = 0, b = 0, l = 0, unit = "cm"),
      rect = element_rect(linewidth = 0.5)
    )

}

########## Violin + Box ##########
VlnBox <- function(data, x, y, cols = NULL) {
  
  count_var <- data %>%
    group_by(data[[x]]) %>%
    tally()
  colnames(count_var) <- c(x, "n")
  
  p <- ggplot(data, aes(data[[x]], data[[y]])) +
    geom_violin(scale = "width", trim = F, aes(fill = data[[x]])) +
    geom_boxplot(outlier.size = 0, width = 0.3, fill = "white") +
    theme_bw() +
    scale_y_continuous(limits = c(0,1)) +
    geom_text(
      data = count_var,
      aes(x = count_var[[x]], y = -Inf, 
          label = paste0("n = ", format(n, big.mark = ",", trim = TRUE)), vjust = -1),
      color = "black", size = 3
    ) +
    ylab(y) +
    scale_y_continuous(labels = scales::comma, expand = c(0.08,0)) +
    scale_fill_manual(values = cols) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      legend.position = "none"
    )
  
  return(p)
  
}

########## Half violin + half box ##########

HalfVlnBox <- function(data, x, y, cols) {
  
  count_var <- data %>%
    group_by(data[[x]]) %>%
    tally()
  colnames(count_var) <- c(x, "n")
  
  p <- ggplot() +
    geom_half_violin(
      data = data, aes(data[[x]], data[[y]], fill = data[[x]]),
      side = "l", show.legend = FALSE, trim = FALSE
    ) +
    geom_half_boxplot(
      data = data, aes(data[[x]], data[[y]], fill = data[[x]]),
      side = "r", outlier.color = NA, width = 0.4, show.legend = FALSE
    ) +
    geom_text(
      data = count_var,
      aes(x = count_var[[x]], y = -Inf, 
          label = paste0("n = ", format(n, big.mark = ",", trim = TRUE)), vjust = -1),
      color = "black", size = 3
    ) +
    ylab(y) +
    scale_y_continuous(labels = scales::comma, expand = c(0.08,0)) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.text = element_text(color = "black")
    )
  
  return(p)
  
}

########### Dimensional reduction by multiple variables ##########

MultiVarDimRed <- function(seurat, var, reduction) {
  
  L_plots <- list()
  for (v in var) {
    
    if (v == "SRSF2") {
      p <- DimPlot(
        final_seurat,
        group.by = v,
        na.value = "lightgrey",
        pt.size = 0.1,
        cells.highlight = list(
          "mutated" = colnames(final_seurat)[which(final_seurat$SRSF2 == "mutated")],
          "wild-type" = colnames(final_seurat)[which(final_seurat$SRSF2 == "wild-type")]
        ),
        sizes.highlight = c(0.1,1.5)
      ) +
        scale_color_manual(
          values = c("lightgrey", "steelblue", "darkred"),
          labels = c("NA", "wild-type", "mutated")
        ) +
        theme_bw() +
        coord_equal() +
        theme(
          panel.grid = element_blank(),
          aspect.ratio = 1,
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 10),
          plot.title = element_blank(),
          #plot.title = element_text(face = 'bold', hjust = 0.5, vjust = 2),
          legend.margin = margin(rep(0,4), "cm"),
          legend.justification = "left",
          legend.text = element_text(size = 10),
          plot.margin = unit(c(0.1, 0.1, 0.1 , 0.1), "cm")
        )
    } else {
      
    p <- DimPlot(
      object = seurat, 
      reduction = reduction, 
      group.by = v, 
      cols = colors[[v]], 
      pt.size = 0.1, 
      shuffle = T, seed = 123,
      na.value = "lightgrey"
      ) +
      theme_bw() +
      coord_equal() +
      theme(
        panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        plot.title = element_blank(),
        #plot.title = element_text(face = 'bold', hjust = 0.5, vjust = 2),
        legend.margin = margin(rep(0,4), "cm"),
        legend.justification = "left",
        legend.text = element_text(size = 10),
        plot.margin = unit(c(0.1, 0.1, 0.1 , 0.1), "cm")
      )
    }
    
    L_plots[[v]] <- p
  }
  
  pgrid <- plot_grid(plotlist = L_plots, ncol = 2, align = "hv")
  return(pgrid)
  
}

########## Jitter plot ##########
JitterPlot <- function(data, x, y, cols, bar, stat_pwc, alternative = NULL, symnum.args) {
  
  data2 <- data %>% group_by(data[x]) %>% summarise(mean = mean(.data[[y]]))
  
  p <- ggplot(data, aes(data[[x]], data[[y]], fill = data[[x]])) +
    geom_jitter(
      size = 2.5, width = 0.2, height = 0, 
      shape = 21, stroke = 0.2, color = "black"
    ) +
    geom_point(
      data = data2,
      aes(x = cohort, y = mean),
      size = 10, color = 'black', fill = "black", shape = '_'
    ) +
    theme_bw() +
    scale_y_continuous(expand = c(0.1,0.1)) +
    ylab(y) +
    scale_fill_manual(values = cols) +
    stat_pwc(
      method = stat_pwc,
      method.args = list(alternative = alternative),
      bracket.shorten = 0.3,
      tip.length = 0,
      vjust = 0,
      label.size = 2.8,
      symnum.args = symnum.args,
      hide.ns = F
    ) +
    theme(
      panel.grid = element_blank(),
      aspect.ratio = 1,
      legend.position = "none", 
      axis.title.x = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 1),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(face = "bold", hjust = 0.5, vjust = 2, size = 12),
      plot.margin = unit(c(1,1,0,0), "cm")
    )
  
  return(p)
  
}
