
########## Create project directory tree ##########

InitializeProjectDirs <- function(project_name) {
  
  if (!dir.exists(paste(getwd(), project_name, sep = "/"))) {
    
    dir.create(paste(getwd(), project_name, sep = "/"))
    dir.create(paste(getwd(), project_name, "code", sep = "/"))
    dir.create(paste(getwd(), project_name, "data", sep = "/"))
    dir.create(paste(getwd(), project_name, "results", sep = "/"))
    
  }
  
  dir_list <- list(
    "01_qc_normalization",
    "02_cell_type_annotation",
    "03_dim_reduction_batch_correction",
    "04_clustering",
    "05_cluster_analysis",
    "06_cell_type_analysis",
    "07_cell_type_deconvolution"
  )
  
  for (d in dir_list) {
    
    if (!dir.exists(paste(getwd(), project_name, "results", d, sep = "/"))) {

      dir.create(paste(getwd(), project_name, "results", d, sep = "/"))
      
      if (d != "00_initialize_project") {
        dir.create(paste(getwd(), project_name, "results", d, "plots", sep = "/"))
        dir.create(paste(getwd(), project_name, "results", d, "tables", sep = "/"))        
      }
    }
  }
}

########## Upload count matrices and set metadata ##########
GetListSeuratObjs <- function(sample_sheet) {
  
  col_metadata <- 
    colnames(sample_sheet)[which(!colnames(sample_sheet) %in% c("chromium_id", "seq_date", "path2tenx", "date_birth", "date_diagnosis"))]
  
  L_seurat <- list()
  for (i in unique(sample_sheet$sample_id)) {
    
    # Log message
    message(paste0(format(Sys.time(), '[%Y-%m-%d %H:%M:%S]'), " Processing sample...",i))
    
    # Import cellranger count matrix
    sample_counts <- Read10X(sample_sheet$path2tenx[which(sample_sheet$sample_id == i)])
    
    # Create Seurat object
    seurat <- CreateSeuratObject(
      counts = sample_counts,
      min.cells = 0,
      min.features = 0
    )
    
    # Set metadata
    for (c in col_metadata) {
      seurat[[c]] <- sample_sheet[[c]][which(sample_sheet$sample_id == i)]
    }
    
    # Add percentage of mitochondrial and ribosomial genes
    seurat$percent_MT <- PercentageFeatureSet(seurat, pattern = "^MT-")
    seurat$percent_RB <- PercentageFeatureSet(seurat, pattern = "^RP[SL]")
    
    # Update cell names with sample_id
    seurat <- RenameCells(seurat, new.names = paste(colnames(seurat), seurat$sample_id, sep = "-"))
    
    # If necessary, merge counts from different transcripts of the same gene
    if (any(duplicated(rownames(seurat))) == TRUE) {
      
      transcript_counts <- GetAssayData(seurat, slot = "count")
      gene_names <- rownames(transcript_counts)
      
      ## get names of genes with multiple entries
      duplicated_gene_names <- unique(gene_names[which(duplicated(gene_names))])
      
      ## log message
      message(
        glue::glue(
          "{format(Sys.time(), '[%Y-%m-%d %H:%M:%S]')} Summing up counts for ",
          "{length(duplicated_gene_names)} genes with multiple entries."
        )
      )
      
      ## go through genes with multiple entries
      L_duplicated <- list()
      for (i in duplicated_gene_names) {
        
        ## extract transcripts counts for current gene
        tmp_counts <- transcript_counts[which(gene_names == i),]
        
        ## make sure at least 2 rows are present
        if ( nrow(tmp_counts) >= 2 ) {
          
          ## sum up counts in each cell
          tmp_counts_sum <- Matrix::colSums(tmp_counts)
          L_duplicated[[i]] <- tmp_counts_sum
        }
      }
      summed_counts <- do.call(rbind, L_duplicated)
      
      ## remove original counts from transcript matrix and the gene name from
      ## the list of gene names
      transcript_counts <- transcript_counts[-c(which(gene_names %in% duplicated_gene_names)),]
      gene_names <- gene_names[-c(which(gene_names %in% duplicated_gene_names))]
      
      ## add summed counts to end of transcript count matrix and current gene
      ## name to end of gene names
      transcript_counts <- rbind(transcript_counts, summed_counts)
      
      # Assign polished counts to Seurat
      seurat@assays$RNA@counts <- transcript_counts
    }
    
    L_seurat[[i]] <- seurat
    
  }
  
  return(L_seurat)

}

########## Check if genes have the same order across seurat objects ##########
sameGeneNames <- function(L_seurat) {
  
  ## create place holder for gene names from each count matrix
  gene_names <- list()
  
  ## go through count matrices
  for ( i in names(L_seurat) ) {
    
    ## extract gene names
    gene_names[[i]] <- rownames(L_seurat[[i]])
  }
  
  ## check if all gene names are the same (using the first matrix as a reference
  ## for the others)
  return(all(sapply(gene_names, FUN = identical, gene_names[[1]])))
}

# IN CASE OF NON IDENTICAL GENES

#for ( i in samples ) {
#  gene_names <- rownames(transcripts$raw[[i]])
#  transcripts$raw[[i]] <- as.data.frame(as.matrix(transcripts$raw[[i]]))
#  transcripts$raw[[i]] <- transcripts$raw[[i]] %>%
#    mutate(gene = gene_names) %>%
#    select(gene, everything())
#} # most likely results in memory issues
#
#transcripts$raw$merged <- 
#  full_join(transcripts$raw$A, transcripts$raw$B, by = "gene") %>%
#  full_join(., transcripts$raw$E, by = "gene") 
#  
#transcripts$raw$merged[ is.na(transcripts$raw$merged) ] <- 0

# check for which samples genes are not identical
#all_gene_names <- c()
#for ( i in samples ) {
#  gene_names <- rownames(transcripts$raw[[i]])
#  all_gene_names[[i]] <- gene_names
#}
#
#identical(names(all_gene_names), samples) # TRUE
#
#for (i in all_gene_names) {
#  print( length(i))
#}
#
#identical(all_gene_names$AML4, all_gene_names$AML5) # TRUE
#identical(all_gene_names$AML4, all_gene_names$sAML1) # FALSE
#identical(all_gene_names$AML4, all_gene_names$AML1) # TRUE
#identical(all_gene_names$AML4, all_gene_names$AML2) # TRUE
#identical(all_gene_names$AML4, all_gene_names$AML3) # TRUE
#identical(all_gene_names$AML4, all_gene_names$hBM1) # TRUE
#identical(all_gene_names$AML4, all_gene_names$hBM2) # TRUE
#identical(all_gene_names$AML4, all_gene_names$hBM3) # TRUE
#
#add_names <- rownames(transcripts$raw$AML4)[which(!rownames(transcripts$raw$AML4) %in% rownames(transcripts$raw$sAML1))] 
#length(add_names) #4683
#
#add_ <- matrix(data = NA, 
#                nrow = length(add_names), 
#                ncol = ncol(transcripts$raw$AML4), dimnames = NULL)
#rownames(add_) <- add_names
#transcripts$raw$AML4 <- rbind(transcripts$raw$AML4, add_)
#dim(transcripts$raw$AML4)