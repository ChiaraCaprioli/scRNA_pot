
########## Set project paths ##########

SetPaths <- function(project_name) {
  path_main <- paste(getwd(), project_name, sep = "/")
  Sys.setenv(PATH_MAIN = path_main)
  
  path_data <- paste0(path_main, "/data")
  Sys.setenv(PATH_DATA = path_data)
  
  path_results <- paste0(path_main, "/results")
  Sys.setenv(PATH_RESULTS = path_results)
  
}

########## Set logger ##########
SetLogger <- function(logfile) {
  
  logfile = logfile
  file_appender = file_appender(
    logfile, 
    append = TRUE, 
    layout = default_log_layout()
  )
  my_logger <- logger(threshold = "INFO", appenders = file_appender)
  return(my_logger)
  
}

########## Update sample sheet with patient age and sequencing batch ##########

UpdateSampleSheet <- function(sample_sheet) {
  
  if (file.exists(paste0(PATH_DATA, "/sample_sheet.tsv"))) {
    
    sample_sheet <- read_tsv(paste0(PATH_DATA, "/sample_sheet.tsv"))
    
    ## patient age
    L_age <- list()
    for (i in unique(sample_sheet$sample_id)) {
      
      x <- as.Date(sample_sheet$date_diagnosis[which(sample_sheet$sample_id == i)], "%d/%m/%Y")
      y <- as.Date(sample_sheet$date_birth[which(sample_sheet$sample_id == i)], "%d/%m/%Y")
      age <- round(as.numeric(x - y) / 365, digits = 2)
      L_age[[i]] <- age
      
    }
    
    df_age <- do.call(rbind, L_age) 
    colnames(df_age) <- "age"
    df_age <- as.data.frame(df_age) %>%
      rownames_to_column(var = "sample_id")
    
    ## sequencing run
    L_seq_batch <- list()
    for (i in unique(sample_sheet$sample_id)) {
      
      x <- str_remove_all(sample_sheet$seq_date[which(sample_sheet$sample_id == i)], pattern = "/")
      seq_batch <- paste0("batch_", x)
      L_seq_batch[[i]] <- seq_batch
      
    }
    
    df_seq_batch <- do.call(rbind, L_seq_batch) 
    colnames(df_seq_batch) <- "seq_batch"
    df_seq_batch <- as.data.frame(df_seq_batch) %>%
      rownames_to_column(var = "sample_id")
    
    ## join and save
    sample_sheet <- plyr::join_all(list(sample_sheet, df_age, df_seq_batch), type = "full")
    
    write_tsv(sample_sheet, paste0(PATH_DATA, "/sample_sheet.tsv"))
    
  } else {
    message("Please upload sample sheet!")
  }
}

########## Table counts according to grouping variables ##########

TableGroupedCounts <- function(data, var1, var2) {
  
  tab <- data %>%
    group_by(data[[var1]], data[[var2]]) %>%
    summarize(count = n()) %>%
    spread(`data[[var2]]`, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    dplyr::select(c(`data[[var1]]`, "total_cell_count", everything())) %>%
    arrange(factor(`data[[var1]]`, levels = levels(data[[var1]])))
  
  colnames(tab)[1] <- var1
  
  return(tab)
  
}

########## Shannon entropy according to grouping variables ##########
ShannonEntropyVar <- function(data, g_var, en_var) {
  data_entropy <- data %>%
    dplyr::select(c(
      unlist(g_var), 
      en_var
    )) %>%
    group_by(
      data[g_var[[1]]],
      data[g_var[[2]]],
      data[en_var]
    ) %>%
    summarise(n = n()) 
  
  df_entropy <- data.frame()
  for ( i in unique(data_entropy[[g_var[[1]]]]) ) {
    g_count <- data_entropy[which(data_entropy[g_var[[1]]] == i),]
    df_entropy <- rbind(
      cbind(
        g1 = unique(g_count[g_var[[1]]]),
        g2 = unique(g_count[g_var[[2]]]),
        entropy = entropy(g_count$n)
      ),
      df_entropy
    )
  }
  df_entropy <- df_entropy %>% 
    arrange(entropy)
  df_entropy[[g_var[[1]]]] <- factor(df_entropy[[g_var[[1]]]], levels = df_entropy[[g_var[[1]]]])
  
  
  p1 <- ggplot(df_entropy, aes(x = df_entropy[[g_var[[1]]]], y = df_entropy[["entropy"]])) +
    geom_point() +
    geom_line(group = g_var[[1]]) +
    labs(y = "Shannon entropy") +
    theme_classic() +
    coord_cartesian(clip = 'off') +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      aspect.ratio = 0.35,
      text = element_text(size = 13, color = 'black'),
      axis.text = element_text(color = 'black'),
      rect = element_rect(linewidth = 0.5),
      axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0))
    )
  
  tab_group_count <- TableGroupedCounts(data = data, var1 = g_var[[1]], var2 = en_var)
  tab_group_count[[g_var[[1]]]] <- factor(tab_group_count[[g_var[[1]]]], levels = df_entropy[[g_var[[1]]]])
  
  p2 <- BarGroupProp(tab_group_count, id.var = g_var[[1]], col = colors[[en_var]]) +
    theme(
      plot.margin = margin(t = 0.3, r = 0, b = 0, l = 0, unit = "cm"),
      axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
    )
  
  plot_save <- p1 / p2
  L_save <- list("df_entropy" = df_entropy, "plot_save" = plot_save)
  return(L_save)
  
}

########## Create list of contrasts by group of interest ##########
GetContrasts <- function(group, levels, sep) {
  
  group <- levels(group)
  cb_list  <- as.list(combn(group, 2, FUN = function(x){paste0(x[1], sep, x[2])}))
  return(cb_list)
  
}

########## Kolmogorov-Smirnov test by multiple contrasts ##########
KsMultipleContrasts <- function(data, dist, var, contrasts) {
  
  # Prepare data
  x <- data %>% select(dist, var)
  x <- x[!is.na(x[[dist]]),]
  x[[dist]] <- x[[dist]]/max(x[[dist]]) #scale 0:1
  
  # Test significance of differences between distributions by Kolmogorov-Smirnov test
  df_ks <- data.frame()
  for (c in contrasts) {
    
    ks <- ks.test(
      x[[dist]][x[[var]] == str_split_1(c, pattern = "_")[1]],
      x[[dist]][x[[var]] == str_split_1(c, pattern = "_")[2]]
    ) # two-sided, pairwise
    
    df_ks <- rbind(
      cbind(
        contrast = c,
        test = "Two-sided Kolmogorov-Smirnov",
        statistic = ks$statistic,
        p_value = ifelse(ks$p.value < 2.2e-16, "< 2.2e-16", round(ks$p.value, digits = 3)) 
      ),
      df_ks
    )
  }
  df_ks <- map_df(df_ks, rev)
  
  # Prepare labels
  labels <- c()
  for (i in 1:nrow(df_ks)) {
    label = paste0(
      df_ks$contrast[[i]], 
      ": D = ", round(as.numeric(df_ks$statistic[[i]]), digits = 2), 
      ", p ", df_ks$p_value[[i]]
    )
    labels[[i]] <- label
  }
  
  # Plot
  p <- ggplot(x, aes(x[[1]], x[[2]], fill = x[[2]])) +
    ggridges::geom_density_ridges(color = "black", size = 0.2) +
    xlab(dist) +
    theme_bw() +
    xlim(0,1) +
    expand_limits(y = c(1, 6)) +
    scale_fill_manual(values = colors[[var]]) +
    annotate(
      geom = "text",
      x = 0.5, y = 5,
      label = str_flatten(labels, collapse = "\n"),
      size = 3
    ) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text = element_text(color = 'black', size = 12),
      axis.title = element_text(color = 'black', size = 12),
      axis.title.y = element_blank(),
      aspect.ratio = 1
    ) 
  
  L <- list("df_ks" = df_ks, "p" = p)
  return(L)
  
}

