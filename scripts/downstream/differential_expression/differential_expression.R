########## Create contrast matrix by group of interest ##########
GetContrastMatrix <- function(group, levels, delim){
  
  # ensure that group levels are unique
  group <- unique(as.character(group))
  
  # make all combinations
  cb  <- combn(group, 2, FUN = function(x){paste0(x[1], "-", x[2])})
  
  # make contrasts
  contrasts <- limma::makeContrasts(contrasts = cb, levels = levels)
  colnames(contrasts) <- gsub("-", delim, colnames(contrasts))
  
  return(contrasts)
  
}