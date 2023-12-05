
########## Barplot proportions of grouped variables ##########

# Requires data = table whose columns contain main grouping variable, total cell count and 
# counts by other grouping variables (order is mandatory)

PlotGroupProp <- function(data, id.var, col = NULL) {

  data_prop <- data %>%
    dplyr::select(-c('total_cell_count')) %>%
    reshape2::melt(id.vars = id.var) 
  data_prop %>% 
    ggplot(aes(data_prop[[id.var]], value)) +
    geom_bar(aes(fill = variable), color = "black", position = 'fill', stat = 'identity', width = 0.7) +
    geom_text(
      data = data,
      aes(x = data[[1]], y = Inf, 
          label = paste0('n = ', format(total_cell_count, big.mark = ',', trim = TRUE)), 
          vjust = 0.5, hjust = -0.2, angle = 90),
      color = 'black', size = 3
    ) +
    scale_y_continuous(name = '% cells', labels = scales::percent_format(), expand = c(0.01,0)) +
    scale_fill_manual(values = c("lightgrey", "#737373")) +
    coord_cartesian(clip = 'off') +
    theme_bw() +
    theme(
      text = element_text(size = 13, color = 'black'),
      aspect.ratio = 1,
      axis.text = element_text(color = 'black'),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      plot.margin = margin(t = 60, r = 0, b = 0, l = 0, unit = 'pt')
    )

}

########## Violin + Box ##########
VlnBox <- function(data, x, y, cols) {
  
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


########### Dimensional reduction ##########
#
#
#function(obj, reduction, group.by, cols, pt.size, shuffle, seed) {
#  
#}
#
#DimPlot(harmony_seurat, reduction = 'UMAP', 
#        group.by = 'lineage', cols = colors$lineage, pt.size = 0.1, shuffle = T, seed = 123) +
#  theme_bw() +
#  coord_equal() +
#  ggtitle("harmony") +
#  theme(panel.grid = element_blank(),
#        aspect.ratio = 1,
#        axis.title = element_text(size = 13),
#        axis.text = element_text(size = 10),
#        plot.title = element_text(face = 'bold', hjust = 0.5, vjust = 2),
#        legend.position = 'bottom'
#  )#