########## Load gene sets ##########
LoadGeneSets <- function(path, ignore_cols = 1){
  x <- scan(path, what="", sep="\n")
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  for(i in 1:ignore_cols){
    y <- lapply(y, `[`, -1) 
  }
  return(y)
}

########## Check gene sets ##########
CheckGeneSetExpr <- function(geneset, expr_matrix, fraction_expressed) {
  
  df <- data.frame()
  for (g in names(geneset)) {
    x = length(which(geneset[[g]] %in% rownames(exprMatrix)))
    y = length(geneset[[g]])
    df <- rbind(
      data.frame(
        signature = g,
        n_gene_set = x,
        fraction_expressed = x / y
      ),
      df
    )
  }
  
  df$keep <- ifelse(df$fraction_expressed >= fraction_expressed, TRUE, FALSE)
  return(df)
  
}

########## Select AUC threshold ##########
# See https://bioconductor.org/packages/devel/bioc/vignettes/AUCell/inst/doc/AUCell.html.

GetMaxAUC <- function(cells_rankings, range, gene_count) {
  L <- list()
  for (i in range){
    aucMaxRank <- ceiling(nrow(cells_rankings)*i)
    df <- data.frame(
      aucMaxRank = aucMaxRank,
      suitable_count = gene_count[["50%"]] >= aucMaxRank
    )
    L[[as.character(i)]] <- df
  }
  x <- do.call(rbind, L)
  aucMaxRank <- max(x$aucMaxRank[which(x$suitable_count == TRUE)])
  message("Using ", rownames(x)[which(x$aucMaxRank == aucMaxRank)], " to compute aucMaxRank")
  return(aucMaxRank)
}

########## Mean AUC by group(s) of interest ##########
GeatMeanAUC <- function(data, group_var, sig) {
  
  # Prepare data
  data <- data %>% dplyr::select(group_var, sig)
  
  if (length(group_var) == 0 | length(group_var) > 2) {
    stop("Check grouping var(s)!")
  }
  
  # Get mean auc by grouping var(s)
  L <- list()
  for (i in setdiff(colnames(data), group_var)) {
    
    if (length(group_var) == 1) {
      x <- aggregate(data[[i]] ~ data[[group_var]], data, mean)
    }
    if (length(group_var) == 2) {
      x <- aggregate(data[[i]] ~ data[[group_var[1]]] + data[[group_var[2]]], data = data, mean)
    } 
    
    L[[i]] <- x
    
  }
  
  sig_score <- do.call(rbind, L) %>% rownames_to_column(var = "signature")
  colnames(sig_score)[2:ncol(sig_score)] <- c(group_var, "score")
  sig_score$signature <- str_split_i(sig_score$signature, pattern = "[.]", 1)
  
  return(sig_score)
  
}