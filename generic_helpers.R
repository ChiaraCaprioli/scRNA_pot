
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
  
  # prepare data
  x <- data_entropy %>%
    group_by(
      data_entropy[g_var[[1]]],
      data_entropy[g_var[[2]]],
      data_entropy[en_var]
    ) %>%
    summarise(n = n()) 
  
  # get entropy by grouping var
  df_entropy <- data.frame()
  for ( i in unique(x[[g_var[[1]]]]) ) {
    g_count <- x[which(x[g_var[[1]]] == i),]
    df_entropy <- rbind(
      cbind(
        g1 = unique(g_count[g_var[[1]]]),
        g2 = unique(g_count[g_var[[2]]]),
        entropy = entropy(g_count$n)
      ),
      df_entropy
    )
  }
  
  df_entropy <- df_entropy %>% map_df(rev)
  return(df_entropy)
  
}
