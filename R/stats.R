

#' Calculates the quantile for a given channel (column).
#' 
#' Compatible with GatingSet, GatingHierarchy, flowSet, flowFrame, cytoset, and cytoframe objects.
#' If a GatingSet, flowSet or cytoset is provided containing multiple samples, all data is 
#' concatenated and the quantile calculated using all of the data.
#' 
#' @param gs GatingSet, GatingHierarchy, flowSet, flowFrame, cytoset or cytoframe object
#' @param channels Character vector specifying the channels (column names) in gs to calculate the quantile for (default: all channels)
#' @param prob A numeric probability value between 0 and 1 or a vector of probabilities (by default, 0.99)
#' @param parentId Optional name of gate to use to filter data before calculating quantile. Only used when a GatingSet or GatingHierarchy is provided
#' @param ... Additional arguments passed to [stats::quantile] 
#' 
#' @returns A named list of quantile values
#' 
#' @examples
#' \dontrun{
#' # Calculate quantile values for the CD4 channels using a flowSet, flowFrame, cytoset, or cytoframe
#' q <- calculate_quantile(fs, "CD4", prob = 0.99)
#' 
#' # Calculate quantile values for the CD4 channels using a GatingSet or GatingHierarchy by filtering for "Live" cells first
#' q <- calculate_quantile(gs, "CD4", prob = 0.99, parentId = "Live")
#' }
#' @export
calculate_quantile <- function(gs, 
                               channels = NULL,
                               probs = 0.99, 
                               parentId = "root",
                               ...){
  
  # Validate inputs
  if (!inherits(gs, c("GatingSet","GatingHierarchy","flowSet","flowFrame","cytoset","cytoframe"))) {
    stop("Argument 'gs' must be a GatingSet, GatingHierarchy, flowSet, flowFrame, cytoset, or cytoframe object")
    
  }else if (inherits(gs, c("GatingSet","GatingHierarchy"))){
    # Get data from GatingSet after filtering by parentId
    # Returns a cytoset
    gs <- gs_pop_get_data(gs, parentId)
  }
  
  # Get expression matrix or matrices
  if (inherits(gs, c("flowSet","cytoset"))){
    # Get list of expression matrices
    expr_mat <- fsApply(gs, function(ff) exprs(ff), simplify = FALSE)
    # Combine by row into one matrix
    expr_mat <- do.call(rbind, expr_mat)
  }else{
    expr_mat <- exprs(gs)
  }
  # Calculate quantile for given channels
  if (!is.null(channels)){
    # Validate channel name exist
    if (!all(channels %in% colnames(gs))){
      stop("Some 'channels' not in gs column names")
    } 
    quantiles <- lapply(channels, function(channel) stats::quantile(expr_mat[,channel], probs=probs, names=FALSE, ...))
    names(quantiles) <- channels
    # Else for all channels
  }else{
    quantiles <- lapply(colnames(gs), function(channel) stats::quantile(expr_mat[,channel], probs=probs, names=FALSE, ...))
    names(quantiles) <- colnames(gs)
  }
  return(quantiles)
}

#' Calculates the range (min, max) for a given channel (column).
#' 
#' Compatible with GatingSet, GatingHierarchy, flowSet, flowFrame, cytoset, and cytoframe objects.
#' If a GatingSet, flowSet or cytoset is provided containing multiple samples, all data is 
#' concatenated and the quantile calculated using all of the data.
#' 
#' @param gs GatingSet, GatingHierarchy, flowSet, flowFrame, cytoset or cytoframe object
#' @param channels Character vector specifying the channels (column names) in gs to calculate the range for (default: all channels)
#' @param parentId Optional name of gate to use to filter data before calculating quantile. Only used when a GatingSet or GatingHierarchy is provided
#' 
#' @returns A named list of range values
#' 
#' @examples
#' \dontrun{
#' # Calculate range values for the CD4 channels using a flowSet, flowFrame, cytoset, or cytoframe
#' r <- calculate_range(fs, "CD4")
#' 
#' # Calculate range values for the CD4 channels using a GatingSet or GatingHierarchy by filtering for "Live" cells first
#' r <- calculate_range(gs, "CD4", parentId = "Live")
#' }
#' @export
calculate_range <- function(gs, 
                            channels = NULL,
                            parentId = "root"){
  
  # Validate inputs
  if (!inherits(gs, c("GatingSet","GatingHierarchy","flowSet","flowFrame","cytoset","cytoframe"))) {
    stop("Argument 'gs' must be a GatingSet, GatingHierarchy, flowSet, flowFrame, cytoset, or cytoframe object")
    
  }else if (inherits(gs, c("GatingSet","GatingHierarchy"))){
    # Get data from GatingSet after filtering by parentId
    # Returns a cytoset
    gs <- gs_pop_get_data(gs, parentId)
  }
  
  # Get expression matrix or matrices
  if (inherits(gs, c("flowSet","cytoset"))){
    # Get list of expression matrices
    expr_mat <- fsApply(gs, function(ff) exprs(ff), simplify = FALSE)
    # Combine by row into one matrix
    expr_mat <- do.call(rbind, expr_mat)
  }else{
    expr_mat <- exprs(gs)
  }
  # Calculate quantile for given channels
  if (!is.null(channels)){
    # Validate channel name exist
    if (!all(channels %in% colnames(gs))){
      stop("Some 'channels' not in gs column names")
    } 
    ranges <- lapply(channels, function(channel) range(expr_mat[,channel]))
    names(ranges) <- channels
    # Else for all channels
  }else{
    ranges <- lapply(colnames(gs), function(channel) range(expr_mat[,channel]))
    names(ranges) <- colnames(gs)
  }
  return(ranges)
}