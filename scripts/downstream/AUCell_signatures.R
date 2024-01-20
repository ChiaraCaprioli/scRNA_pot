########## Create gene sets from signatures ##########

##' Create gene sets from signatures and save as .gmt files.

##' @title CreateGeneSets
##' @param sig_sheet A sheet with signatures of interest. This needs the following fields: 
##' "signature" (mandatory, i.e. the name of the signature to include); 
##' "genes" (mandatory, i.e. gene symbol included in the signature when manually specified);
##' "notes" (mandatory, i.e. the source for the signature, including MSigDB);
##' "collection" (optional, i.e. a name for a collection that groups different signatures).

CreateGeneSets <- function(sig_sheet, path_save) {
  
  # Get signature sheet
  sig <- read_delim(sig_sheet, show_col_types = FALSE)
  
  # Get signatures from MSigDB if necessary
  if (any(sig$notes == "MSigDB")) {
    
    message("Getting signatures from MSigDB...")
    msig_db <- msigdbr(species = "Homo sapiens")
    
  }
  
  message("Creating gene set of signatures and/or collections...")
  
  if ("collection" %in% colnames(sig)) {
    
    collection <- unique(sig$collection)
    
    # create GeneSet and save as .gmt
    for (c in collection) {
      
      sig_collection <- sig$signature[which(sig$collection == c)]
      L_sig <- list()
      
      for (s in sig_collection) {
        
        if (any(msig_db$gs_name %in% s)) {
          geneIds <- unique(msig_db$gene_symbol[which(msig_db$gs_name == s)])
        } else {
          geneIds <- str_split_1(sig$genes[which(sig$signature == s)], pattern = ",")
        }
        
        gs <- GeneSet(geneIds = geneIds)
        gs@setName <- s
        L_sig[[s]] <- gs
        
      }
      
      # create collection
      gs_collection <- GeneSetCollection(L_sig)
      
      # save as .gmt
      toGmt(gs_collection, paste0(path_save, "/", c, ".gmt"))
      
    }
  }
  
  if (!"collection" %in% colnames(sig)) {
    
    signature <- unique(sig$signature)
    
    # create GeneSet and save as .gmt
    for (s in signature) {
      
      geneIds <- str_split_1(sig$genes[which(sig$signature == s)], pattern = ",")
      gs <- GeneSet(geneIds = geneIds)
      gs@setName <- s
      toGmt(gs, paste0(path_save, "/", s, ".gmt"))
      
    }
  }
  
  message("Done!")
  
}

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
GetMeanAUC <- function(data, group_var, sig) {
  
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

########## Z-score of AUC by group(s) of interest ##########
ZGroupAUC <- function(data, group_var, sig) {
  
  # Prepare data
  data <- data %>% dplyr::select(group_var, sig)
  
  L <- list()
  for (i in setdiff(colnames(data), group_var)) {
    
    m_group <- aggregate(data[[i]] ~ data[[group_var]], data, mean)
    m_all <- mean(data[[i]])
    sd_all <- sd(data[[i]])
    z <- (m_group[[2]] - m_all) / sd_all
    
    df_z <- data.frame(
      group_var = m_group[[1]],
      z_score = z
    )
    
    L[[i]] <- df_z
    
  }
  
  sig_zscore <- do.call(rbind, L) %>% rownames_to_column(var = "signature")
  colnames(sig_zscore)[2:ncol(sig_zscore)] <- c(group_var, "z_score")
  sig_zscore$signature <- str_split_i(sig_zscore$signature, pattern = "[.]", 1)
  
  return(sig_zscore)
  
}