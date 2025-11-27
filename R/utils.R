
#' Extract experiment metadata from .fcs filenames and add it to phenoData.
#' 
#' This function uses regular expressions to find patterns in the filenames.
#' 
#' @param fs A flowSet or GatingSet object
#' @param ... Regular expression patterns to extract. Provide as id="regex pattern" where "id" will be added as a column in phenoData with the returned metadata as values in that column
#' 
#' @returns The flowSet or GatingSet with the metadata added. Can view added metadata by calling pData(fs).
#' 
#' @examples
#' \dontrun{
#' # Add a group column identifying samples and controls 
#' fs <- read_filenames_to_pdata(fs, group = "(?i)scc|fmo|unstained|sample")                              
#' }
#' @export
read_filenames_to_pdata <- function(fs, ...){
  search_terms <- list(...)
  # Validate input
  if (length(search_terms) < 1){
    stop("Provide some search keywords")
  }
  if (any(names(search_terms) == "") | any(is.null(names(search_terms)))){
    stop("Invalid search keywords")
  }
  # Get pData, will be used multiple times
  pd <- pData(fs)
  # Search for each key in pData
  for (key in names(search_terms)){
   res <- lapply(pd, function(x) stringr::str_extract(x, pattern=search_terms[[key]]))
   
   if (any(is.na(res$name))) warning(sprintf("Could not find a match for %s=%s in all filename.",key, search_terms[[key]]))
   # Add to pData
   pData(fs)[[key]] <- res$name
  }
  return(fs)
}

#' Searches the phenoData for provided keywords and returns the matching samples.
#' If a match is found, returns the associated flowFrame when a flowSet is provided, or a GatingHierarchy when a GatingSet is provided.
#' If multiple matches are found, a flowSet/GatingSet is returned containing those matches. If no matches, the original data is returned.
#' 
#' Provide keywords as col="keyword", where col is the name of a column in phenoData. 
#' Multiple keywords can be provided. When match_all = TRUE, a sample must match all the keywords to be returned.
#' When FALSE, samples that match any keyword are returned.
#' 
#' @param fs A flowSet or GatingSet
#' @param ... Keywords to search for. Should be provided as named arguments where the name is the column to search in the phenoData and the value is the search keywords  
#' @param use_regex Whether to search using regular expressions (default: TRUE)
#' @param ignore_case Whether the search is case sensitive. Only applies when use_regex = TRUE (default: TRUE)
#' @param match_all Whether to return samples that match all or any keywords (default: TRUE, must match all)
#' @param verbose Whether to print the samples that match
#' 
#' @returns A flowFrame, flowSet, GatingSet, or GatingHierarchy depending on what was provided
#' 
##' @examples
#' \dontrun{
#' # Get all samples in group "FMO"
#' fs_subset <- get_samples_by_keyword(fs, group = "FMO")
#' }
#' 
#' @export
get_samples_by_keyword <- function(fs,
                                    ...,
                                    use_regex=TRUE,
                                    ignore_case=TRUE,
                                    match_all=TRUE,
                                    verbose=TRUE){
  search_terms <- list(...)
  # Validate input
  if (length(search_terms) < 1){
    stop("Provide some search keywords")
  }
  if (any(names(search_terms) == "") | any(is.null(names(search_terms)))){
    stop("Invalid search keywords")
  }
  # Get pData to search
  pd <- pData(fs)
  
  if (use_regex){
    sr <- lapply(names(search_terms), function(name){
      tryCatch(grepl(search_terms[[name]], pd[[name]], ignore.case=ignore_case), error = function(e) FALSE)
      })
  }else{
    sr <- lapply(names(search_terms), function(name){
      tryCatch(pd[[name]]==search_terms[[name]], error = function(e) FALSE)
      })
  }
  if (match_all){
    # Logical "|" to include those that match all keywords
    search_res <- Reduce("&", sr)
  }else{
    # Logical "&" to include those that match any keywords
    search_res <- Reduce("|", sr)
  }
  # Count number of rows meeting search
  count_res <- sum(search_res, na.rm = TRUE)
  # Get indices of matching
  idx <- which(search_res)
  # If one result found
  if (count_res==1){
    ff <- fs[[idx,]]
    if (verbose) message(cat(count_res," samples found:\n",paste0(pd[idx,"name"], collapse="\n"),sep=""), appendLF=TRUE)
  }else if (count_res > 1){
    ff <- fs[idx,]
    if (verbose) message(cat(count_res," samples found:\n",paste0(pd[idx,"name"], collapse="\n"),sep=""), appendLF=TRUE)
  }else{
    ff <- fs
    if (verbose) message(cat("No samples found",sep=""), appendLF=TRUE)
  }
  return(ff)  
}


#' Returns the spill matrix contained in a flowFrame. If none is present, returns NULL.
#' 
#' @param ff A flowFrame object
#' 
#' @returns The spill matrix for the flowFrame, if found, or else NULL
get_spill_matrix <- function(ff){
  # Get spill matrix if present - returns a named list
  # Equivalent to keyword(x, c("spillover", "SPILL", "$SPILLOVER")
  spill <- spillover(ff)
  # Remove NULL values
  spill <- spill[!sapply(spill,is.null)]
  # Return spill if found, else NULL
  if (length(spill)==1){
    return (spill[[1]])
  }else{
    message(sprintf("Could not identify spill matrix for: %s",identifier(ff)))
    return (invisible(NULL))
  }
}

#' Apply compensation from spill matrices in provided in .fcs files.
#' Compensation should be applied before transforming data and before renaming or filtering channels.
#' 
#' @param fs A flowSet object
#' 
#' @returns The flowSet with compensation applied.
#' 
#' @export
apply_compensation_from_fcs <- function(fs){
  # Get named list of spillover from each flowframe
  comp <- fsApply(fs, get_spill_matrix, simplify=FALSE)
  # Apply to flowSet
  fs_comp <- compensate(fs, comp)
  return(fs)
}

#' Returns the marker names contained in flowSet.
#' The marker names are typically user-defined when the .fcs files are created and can include descriptive names like "CD45 Y585-PE-A".
#' You can call markernames(fs) to view the marker names contained in the flowSet.
#' Does not alter scatter (FCS-A, SSC-A, etc.) or time channels.
#' 
#' Optionally, regular expression(s) can be provided to extract substrings from marker names.
#' When multiple patterns are provided, the substrings will be concatenated together in the order provided.
#' This can be used to avoid duplicate channel names. For example, when both area and height are recorded for a given marker, it may be useful to distinguish them using the pattern "-A|-H".
#' If the provided pattern does not match, the full marker name is returned unless partial matches are permitted.
#' 
#' If desired, a default name can be provided for marker names without a match. This can be used to label unused channels
#' 
#' @param fs A flowSet containing multiple .fcs files
#' @param ... Any number of regular expression patterns (strings) to extract from marker names.
#' @param sep The separator to use when concatenating strings (default: "")
#' @param default_name Optional name for marker names where no patterns match (default: none)
#' @param allow_partial_match Whether to use marker names when not all regular expression patterns match (default: TRUE)
#' @param verbose Whether to print updated markers (default: TRUE)
#' 
#' @returns A list of marker names for their corresponding channels. This can be applied to fs
#'  using [set_channels_from_markers()].
#'  
#' @examples
#' \dontrun{
#' # Get marker names
#' markers <- get_marker_names(fs)
#' }
#' @export
get_marker_names <- function(fs, ..., sep = "", default_name = NULL, allow_partial_match = TRUE, verbose = TRUE){
  patterns <- list(...)

  # Get marker names from flowSet
  markers <- markernames(fs)
  # Use patterns to find substrings in markers
  if (length(patterns) > 0){
    for (name in names(markers)){
      marker <- lapply(patterns, function(pattern) stringr::str_extract(markers[[name]], pattern=pattern))
      # If all patterns match
      if (!any(is.na(marker))){
        # Update name
        markers[[name]] <- paste0(marker, collapse=sep)
      # Else if some patterns match and partial matches allowed
      }else if (!all(is.na(marker)) & allow_partial_match){
        # Remove NA before updating name
        marker <- marker[!is.na(marker)]
        markers[[name]] <- paste0(marker, collapse=sep)
      # Else if default names provided
      }else if (!is.null(default_name)){
        markers[[name]] <- default_name
      }
      # Else marker name used as is
      
      if (verbose) message(sprintf("%s -> %s", name, markers[[name]]))
    }
  }else if (verbose){
    message("No regular expression patterns provided. Returning marker names as provided.")
  }
  return (markers)
}

#' Set channels (column names) from the marker names in a flowSet. The original channel names
#' become the marker names. The marker names can be provided as a list, such as what is 
#' returned from [get_marker_names()]. If not provided, the marker names contained in the flowSet are used.
#' 
#' @param fs A flowSet object
#' @param markers A list of marker names to set as channel names
#' 
#' @returns The flowSet with updated channel and marker names
#' @export
set_channels_from_markers <- function(fs, markers = NULL){
  # Get marker names from flowSet if not provided
  if (is.null(markers)){
    markers <- markernames(fs)
  }
  # Update column names
  colnames(fs) <- lapply(colnames(fs), function(name) tryCatch(markers[[name]], error=function(e) name))
  # Use old column names as markernames
  markernames(fs) <- setNames(names(markers),markers)
  
  return(fs)
}

#' Filter channels (column names) in flowSet using regular expressions.
#' 
#' @param fs The flowSet to filter.
#' @param keep A regular expression used to select channels to keep in the flowSet.
#' @param exclude A regular expression used to select channels to remove from the flowSet.
#' @param verbose Whether to print the remaining channels (default: TRUE)
#' @returns The filtered flowSet. 
#' 
#' @examples
#' \dontrun{
#' # Keep only scatter and time channels
#' fs <- filter_channels(fs, keep="FSC|SSC|Time")
#' }
#' 
#' @export
filter_channels <- function(fs, keep=NULL, exclude=NULL, verbose=TRUE){
  # Get column names
  channels <- colnames(fs)
  
  # If neither keep or exclude provided, throw warning, no filtering
  if (all(is.null(c(keep,exclude)))){
    warning("Must provide a value to one of keep or exclude (nothing filtered).")
  # If keep provided
  }else if (!is.null(keep)){
    channels <- channels[grepl(keep, channels)]
  # If exclude provided
  }else if (!is.null(exclude)){
    channels <- channels[!grepl(exclude, channels)]
  # both provided, throw warning, no filtering
  }else{
    warning("Provide one of keep or exclude, not both (nothing filtered).")
  }
  # Remove from flowSet
  fs <- fs[,channels]
  
  if (verbose){
    # Print channels retained after filtering
    #channels <- paste(channels, collapse="\n")
    cat("Filterning complete","Channels kept:",channels, sep="\n")
  }
  return(fs)
}

#' Estimate reasonable limits for each channel of a gate in a Gatingset.
#'
#' The values are calculated for the entire dataset, allowing for consistent limits for plotting.
#'
#' @param gs A GatingSet object
#' @param limits A numeric vector of the quantiles (min, max) to use to estimate the limits. If not provided, the data range is used instead (default: data range)
#'
#' @return A named list where each name is a gate and each element is the corresponding limits for the channel(s) of that gate
#'   
#' @export
estimate_channel_limits <- function(gs, limits = NULL) {
  # Validate inputs
  if (is.null(limits)){
    probs <- NULL
  }else if (inherits(limits, "numeric")){
    if (length(limits)==2){
      probs <- limits
    }else{
      stop("Argument 'limits' be have length 2")
    }
  }else{
    stop("Argument 'limits' must be a numeric")
  }
  # Validate inputs
  if (!inherits(gs, "GatingSet")) {
    stop("Input 'gs' must be a GatingSet object")
  }
  
  # Get all nodes (gates) in the hierarchy
  nodes <- gs_get_pop_paths(gs, path = 2)
  # Remove root from list
  nodes <- nodes[nodes != "root"]
  # Nodes whose parent is root
  nodes <- strsplit(nodes,"/")
  # Add root as parent to any single node without a parent
  nodes <- lapply(nodes, function(node){if(length(node)==1) c("root",node[[1]]) else node})
  
  # Initialize empty list to store limits for each gate
  node_limits <- vector("list",length(nodes))
  names(node_limits) <- vapply(nodes, `[`, 2, FUN.VALUE=character(1))
  
  for (node in nodes){
    
    # Get gate information
    gate <- gs_pop_get_gate(gs, node[[2]])
    
    # If returned gate is a list of gates, one for each sample, use the first
    if (inherits(gate,"list")){
      gate <- gate[[1]]
    }
    # Determine which parameters/channels are used in this gate
    if (inherits(gate, c("polygonGate","rectangleGate","ellipsoidGate","quadGate"))) {
      params <- parameters(gate)
    } else if (inherits(gate, c("booleanFilter","logical"))) {
      # Skip boolean gates as they don't have associated channels
      next
    } else {
      # For other gate types, try to extract parameters
      params <- tryCatch(parameters(gate), error = function(e) NULL)
      if (is.null(params) | length(params) == 0) {
        next
      }
    }
    # Get limits as range or quantile
    if (is.null(probs)){
      lim <- calculate_range(gs, channels = c(params[[1]], params[[2]]), parentId = node[[1]])
    }else{
      lim <- calculate_quantile(gs, channels = c(params[[1]], params[[2]]), probs = probs, parentId = node[[1]])
    }
    # Store values
    node_limits[[node[2]]] <- lim
  }
  return (node_limits)
}