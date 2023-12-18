
########## Compute clusters by Leiden algorithm ##########

GetLeidenClusters <- function(seurat, assay, k_list, res_list, seed, reduction, ndims, path_log) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  DefaultAssay(seurat) <- "RNA"
  
  res_L <- list() 
  for (k in k_list) {
    
    k_start_time <- Sys.time()
    
    for (r in res_list) {
      
      info(my_logger, paste0("Evaluating k ", k, " and res ", r))
      
      seurat <- FindNeighbors(
        seurat,
        reduction = reduction, 
        k.param = k,
        dims = 1:ndims,
        do.plot = F
        )
      
      seurat <- seurat %>%
        FindClusters(
          resolution = r, 
          algorithm = 4, # 4 = leiden
          method = "igraph", 
          random.seed = seed
          ) 
      
      res_L[[paste0("RNA_snn_k_", k, "_res_", r)]] <- seurat@meta.data[paste0("RNA_snn_res.", r)]
      
    }
    
    k_end_time <- Sys.time()
    
    info(my_logger, 
         paste0(
           "k ", k, " successfully done in ", 
           floor(as.numeric(k_end_time - k_start_time, units = "mins")), " mins"
           )
    )
    
   # write.table(
   #   readLines(logfile), 
   #   file = paste(path_log, logfile, sep = "/"), 
   #   sep = "\t", 
   #   row.names = FALSE, 
   #   col.names = FALSE
   # )
    
    }
  
  cluster_all <- do.call(cbind, res_L) 
  colnames(cluster_all) <- names(res_L)
  return(cluster_all)
  
}
