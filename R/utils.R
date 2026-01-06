
#' Extract metadata from the "names" column in phenoData and add it as a new
#' phenoData column.
#'
#' This function uses regular expressions to match patterns in the filenames,
#' extracts those matches and adds them to a new column in phenoData.
#'
#' @param fs A flowSet or related object with a phenoData slot accessible using `pData(fs)`
#' @param ... Regular expression patterns to extract. Provide as `id="regex pattern"` where "id" will be added as a column in phenoData with the returned metadata as values in that column
#'
#' @returns The object with the metadata added. Can view added metadata using `pData(fs)`.
#'
#' @examples
#' \dontrun{
#' # Add a group column identifying samples and controls
#' fs <- read_names_to_pdata(fs, group = "(?i)scc|fmo|unstained|sample")
#' }
#' @importFrom stringr str_extract
#' @import flowCore
#' @import flowWorkspace
#' @export
read_names_to_pdata <- function(fs, ...){
  # Validate input
  if (!inherits(fs, c("flowSet","cytoset","GatingSet"))){
    stop("Provide a flowSet, cytoset, or GatingSet object")
  }
  search_terms <- list(...)
  # Validate input
  if (length(search_terms) < 1){
    stop("Provide some regular expression patterns to match")
  }
  if (any(names(search_terms) == "") | any(is.null(names(search_terms)))){
    stop("All regular expression patterns must be provided as named arguments")
  }
  # Get pData, will be used multiple times
  pd <- pData(fs)
  # Search for each key in pData
  for (key in names(search_terms)){
   res <- sapply(rownames(pd), function(x) str_extract(x, pattern=search_terms[[key]]))

   if (all(is.na(res))) warning(sprintf("Could not find a match for %s = '%s'",key, search_terms[[key]]))
   # Add new column
   pd[[key]] <- res
  }
  # Update in fs
  pData(fs) <- pd
  return(fs)
}

#' Searches the phenoData for provided keywords and returns the matching samples.
#'
#' Provide keywords as col="keyword", where col is the name of a column in phenoData.
#' Multiple keywords can be provided. When match_all = TRUE, a sample must match all the keywords to be returned.
#' When FALSE, samples that match any keyword are returned.
#'
#' @param fs A flowSet or related object with a phenoData slot accessible using `pData(fs)`
#' @param ... Keywords to search for. Should be provided as named arguments where the name is the column to search in the phenoData and the value is the search keywords
#' @param use_regex Whether to search using regular expressions (default: TRUE)
#' @param ignore_case Whether the search is case sensitive. Only applies when use_regex = TRUE (default: TRUE)
#' @param match_all Whether to return samples that match all or any keywords (default: TRUE, must match all)
#' @param verbose Whether to print the samples that match
#' @returns If multiple matches are found, the returned object is the same class as
#' provided. If a single match, the object contained in the provided object. For example,
#' if a flowSet is provided, a flowFrame is returned.
#' @examples
#' \dontrun{
#' # Get all samples in group "FMO"
#' fs_subset <- get_samples_by_keyword(fs, group = "FMO")
#' }
#' @import flowCore
#' @import flowWorkspace
#' @export
get_samples_by_keyword <- function(fs,
                                    ...,
                                    use_regex=TRUE,
                                    ignore_case=TRUE,
                                    match_all=TRUE,
                                    verbose=TRUE){
  # Validate input
  if (!inherits(fs, c("flowSet","cytoset","GatingSet"))){
    stop("Provide a flowSet, cytoset, or GatingSet object")
  }
  search_terms <- list(...)
  # Validate input
  if (length(search_terms) < 1){
    stop("Provide some search keywords")
  }
  if (any(names(search_terms) == "") | any(is.null(names(search_terms)))){
    stop("Search keywords must be provided as named arguments")
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
    if (verbose) message(cat("No samples found. Original object returned",sep=""), appendLF=TRUE)
  }
  return(ff)
}


#' Returns the spill matrix contained in a flowFrame.
#'
#' If none is present, returns NULL.
#'
#' @param ff A flowFrame or cytoframe object
#' @returns The spill matrix for the flowFrame, if found, or else NULL
#'
#' @import flowCore
#' @import flowWorkspace
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

#' Apply compensation from spill matrices provided in .fcs files.
#' Compensation should be applied before transforming data and before renaming or filtering channels.
#'
#' @param fs A flowSet, flowFrame, cytoset, or cytoframe object
#' @returns The flowSet with compensation applied.
#' @import flowCore
#' @import flowWorkspace
#' @export
apply_compensation_from_fcs <- function(fs){
  # Validate input
  if (!inherits(fs, c("flowSet", 'flowFrame', "cytoset", "cytoframe"))){
    stop("Only flowSet, flowFrame, cytoset, and cytoframe objects are compatible with this function")
  }

  if (inherits(fs, c("flowSet", "cytoset"))){
    # Get named list of spillover from each flowframe
    comp <- fsApply(fs, get_spill_matrix, simplify=FALSE)
  }else{
    comp <- get_spill_matrix(fs)
  }
  # Apply compensation
  fs <- compensate(fs, comp)
  return(fs)
}

#' Returns the marker names stored in .fcs files.
#'
#' The marker names are typically user-defined when the .fcs files are created
#' and can include descriptive names like "CD45 Y585-PE-A". You can use markernames(fs)
#' to view the current marker names. Does not alter scatter (FCS-A, SSC-A, etc.)
#' or time channels.
#'
#' Optionally, regular expression(s) may be used to extract substrings from marker names.
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
#' @returns A list of marker names for their corresponding channels. This can be applied to fs
#'  using [set_channels_from_markers()].
#' @examples
#' \dontrun{
#' # Get marker names
#' markers <- get_marker_names(fs)
#' }
#' @importFrom stringr str_extract
#' @import flowCore
#' @import flowWorkspace
#' @export
get_marker_names <- function(fs,
                             ...,
                             sep = "",
                             default_name = NULL,
                             allow_partial_match = TRUE,
                             verbose = TRUE){
  # Validate input
  if (!inherits(fs,c("flowSet","flowFrame","cytoset","cytoframe","GatingSet","GatingHierarchy"))){
    stop("Only flowSet, flowFrame, cytoset, and cytoframe, GatingSet, or GatingHierarchy objects are compatible with this function")
  }
  patterns <- list(...)

  # Get marker names from flowSet
  markers <- markernames(fs)
  # Use patterns to find substrings in markers
  if (length(patterns) > 0){
    for (name in names(markers)){
      marker <- lapply(patterns, function(pattern) str_extract(markers[[name]], pattern=pattern))
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
    message("No regular expression patterns provided. Returning marker names as is")
  }
  return (markers)
}

#' Set channels from the marker names in a flowSet or GatingSet.
#'
#' The original channel names become the marker names. The marker names can be provided as a list, such as what is
#' returned from [get_marker_names()]. If not provided, the marker names contained in the flowSet are used.
#'
#' @param fs A flowSet or GatingSet object
#' @param markers A list of marker names to set as channel names. If NULL, the stored marker names are used
#' @param set_markers_as_channels Whether to set the old channel names as marker names (default: FALSE)
#' @returns The flowSet with updated channel and marker names
#' @import flowCore
#' @import flowWorkspace
#' @importFrom stats setNames
#' @export
set_channels_from_markers <- function(fs,
                                      markers = NULL,
                                      set_markers_as_channels = FALSE){
  # Validate input
  if (!inherits(fs,c("flowSet","flowFrame","cytoset","cytoframe","GatingSet","GatingHierarchy"))){
    stop("Only flowSet, flowFrame, cytoset, and cytoframe, GatingSet, or GatingHierarchy objects are compatible with this function")
  }
  # Get marker names from flowSet if not provided
  if (is.null(markers)){
    markers <- markernames(fs)
  }
  # Update column names
  cols <- sapply(colnames(fs), function(name) tryCatch(markers[[name]], error=function(e) name))

  # Check for duplicate names (throws error)
  dupl <- cols[duplicated(cols)]
  # Number duplicates 1,2,...
  if (length(dupl) > 0){
    dupl <- sapply(seq_along(dupl), function(i) paste0(dupl[[i]],".",i))
    cols[duplicated(cols)] <- dupl
  }
  # Add updated column names
  colnames(fs) <- cols
  # Make old column name the marker name
  if (set_markers_as_channels){
    markernames(fs) <- setNames(names(markers),markers)
  }

  return(fs)
}

#' Filter channels in set using regular expressions.
#'
#' @param fs The flowSet, flowFrame, cytoset, or cytoframe object
#' @param keep A regular expression used to select channels to keep in the set
#' @param exclude A regular expression used to select channels to remove from the set
#' @param filter_scatter Whether to allow filtering of scatter channels. When TRUE, the FSC/SSC channels can be removed from the set (default: FALSE)
#' @param ignore_case Whether to ignore case in the regular expression such that "fl1-a" matches "FL1-A" etc (default: TRUE)
#' @param verbose Whether to print the remaining channels (default: TRUE)
#' @returns The filtered flowSet.
#' @examples
#' \dontrun{
#' # Keep only scatter and time channels
#' fs <- filter_channels(fs, keep="FSC|SSC|Time")
#' }
#' @import flowCore
#' @import flowWorkspace
#' @export
filter_channels <- function(fs,
                            keep=NULL,
                            exclude=NULL,
                            filter_scatter=FALSE,
                            ignore_case=TRUE,
                            verbose=TRUE){
  # Validate input
  if (!inherits(fs, c("flowSet", 'flowFrame', "cytoset", "cytoframe"))){
    stop("Only flowSet, flowFrame, cytoset, and cytoframe objects are compatible with this function")
  }
  # Get column names
  channels <- colnames(fs)

  # For debugging method imports - may remove in future
  if (is.null(channels)){
    stop("Failed to get channel names")
  }
  # Split scatter and fluorescent channels - only filter fluorescent
  if (!filter_scatter){
    scatter <- channels[grepl("fsc|ssc", channels, ignore.case = TRUE)]
    channels <- channels[!grepl("fsc|ssc", channels, ignore.case = TRUE)]
  }else{
    scatter <- NULL
  }

  # If neither keep or exclude provided, throw warning, no filtering
  if (all(is.null(c(keep,exclude)))){
    warning("Must provide a value to one of keep or exclude (nothing filtered)")
  # If keep provided
  }else if (!is.null(keep)){
    channels <- channels[grepl(keep, channels, ignore.case = ignore_case)]
  # If exclude provided
  }else if (!is.null(exclude)){
    channels <- channels[!grepl(exclude, channels, ignore.case = ignore_case)]
  # both provided, throw warning, no filtering
  }else{
    warning("Provide one of keep or exclude, not both (nothing filtered)")
  }
  # Concat with scatter (or NULL)
  channels <- c(scatter, channels)
  # Remove from flowSet
  fs <- fs[,channels]

  if (verbose){
    # Print channels retained after filtering
    channels <- paste0(channels, collapse="\n")
    message(sprintf("Channels kept:\n%s",channels))
  }
  return(fs)
}

#' Helper function to flatten flowSet into flowFrame
#'
#' @param fs A flowSet
#' @returns A flowFrame
#' @import flowCore
flatten_to_flowframe <- function(fs){
  # Get list of expression matrices
  expr_mat <- fsApply(fs, function(ff) exprs(ff), simplify = FALSE)
  # Combine by row into one matrix
  expr_mat <- do.call(rbind, expr_mat)
  # Make into flowframe and return
  return (flowFrame(expr_mat))
}

#' Helper function to flatten cytoset into cytoframe
#'
#' @param cs A cytoset
#' @returns A cytoframe
#' @import flowWorkspace
flatten_to_cytoframe <- function(cs){
  # Get list of expression matrices
  expr_mat <- fsApply(cs, function(cf) exprs(cf), simplify = FALSE)
  # Combine by row into one matrix
  expr_mat <- do.call(rbind, expr_mat)
  # Make into cytoframe and return
  return (cytoframe(expr_mat))
}

#' Flatten a set by group
#'
#' Groups frames by a column in phenoData then concatenates each group into a single
#' frame. Intended use is when the set includes replicates that you would like to combine.
#'
#' The resultant frames will have their sampleNames attribute as the value in groupby.
#' All phenoData that is common to the group will be retained, and unique values
#' (like filenames) are discarded.
#'
#' In future updates this function will return the same object class (flowSet or
#' cytoset) as provided. Currently, a flowSet is returned.
#'
#' @param fs A flowSet or cytoset object
#' @param groupby A column in phenoData to use to group frames
#' @param rm_na Whether the exclude frames with NA values in groupby column (default: TRUE)
#' @returns A flowSet containing flattened frame(s)
#' @import flowCore
#' @import flowWorkspace
#' @export
flatten_by_group <- function(fs, groupby, rm_na = TRUE){

  # Validate inputs
  if (!inherits(fs, c("flowSet","cytoset"))) {
    stop("Argument 'fs' must be a flowSet or cytoset object")
  }

  pd <- pData(fs)
  if (!(groupby %in% colnames(pd))){
    stop("Argument 'groupby' must be a column in phenoData")
  }
  # Get unique values in groupby column
  groups <- unique(pd[[groupby]])

  # Remove NA
  if (rm_na) {
    groups <- groups[!is.na(groups)]
    # Check that elements remain
    if (length(groups) < 1) stop(sprintf("Only NA found in column %s",groupby))
  }
  message(sprintf("Groups: %s",paste0(groups, collapse=", ")))

  # Make list of flattened flowframes
  fs <- lapply(groups, function(group){
    # Select frames by group
    fs_sub <- fs[pd[[groupby]]==group]
    # Flatten to frame
    flatten_to_flowframe(fs_sub)
  })
  # Make list of frames into flowSet
  fs <- flowSet(fs)
  # Change default sampleNames
  sampleNames(fs) <- groups

  # Make list of pData for new flowSet
  pd <- lapply(groups, function(group){
    # Select pData by group
    pd_sub <- pd[pd[[groupby]]==group,]
    # Only keep rows common to group
    pd_sub <- subset(pd_sub,select = lengths(sapply(pd_sub, unique))==1)
    # Select first row to keep (all others are the same)
    pd_sub[1,]
  })
  # Make into a single dataframe
  pd <- do.call(rbind, pd)
  # Make rownames as groups (must match sampleNames(fs) or else error)
  rownames(pd) <- groups
  # Now set pData
  pData(fs) <- pd

  return (fs)
}

#' Estimate reasonable channel limits for a gate in a Gatingset.
#'
#' The values are calculated from all samples in the GatingSet, filtered by the gate's
#' parent gate. Intended use if for setting consistent limits for all samples when plotting.
#'
#' @param gs A GatingSet object
#' @param gateId Character name of a gate
#' @param limits A numeric vector of the quantiles (min, max) to use to estimate the limits. If not provided, the data range is used instead (default: data range)
#' @return A named list where each name is a gate and each element is the corresponding limits for the channel(s) of that gate
#' @import flowWorkspace
#' @export
calculate_channel_limits <- function(gs, gateId, limits = NULL) {
  # Validate inputs
  if (!inherits(gs, "GatingSet")) {
    stop("Input 'gs' must be a GatingSet object")
  }
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

  # Get gate
  gate <- gs_pop_get_gate(gs, gateId)
  # Get parent name
  parentId <- gs_pop_get_parent(gs, gateId)

  # If returned gate is a list of gates, one for each sample, use the first
  if (inherits(gate,"list")){
    gate <- gate[[1]]
  }
  # Determine which parameters/channels are used in this gate
  # tryCatch to deal with boolean or logical gates that may not have associated channels
  params <- tryCatch(parameters(gate), error = function(e) NULL)

  # Get limits as range or quantile
  # If params is NULL the limits are calculated for all channels
  # This will prevent failure here but the appropriate channel needs to be selected later
  if (is.null(probs)){
    lim <- calculate_range(gs, channels = params, parentId = parentId)
  }else{
    lim <- calculate_quantile(gs, channels = params, probs = probs, parentId = parentId)
  }
  return(lim)
}

#' Randomly sample from frames in a set, returning n events per frame.
#' If the number of events in the frame < n, all events are retained.
#'
#' @param fs Flowset or cytoset object
#' @param n Number of events to sample (default: 25000)
#' @return The downsampled set
#' @import flowCore
#' @import flowWorkspace
#' @export
downsampleFrames <- function(fs, n = 25000){
  # Validate input
  if (!inherits(fs, c("flowSet", "cytoset"))){
    stop("Only flowSet and cytoset objects are compatible with this function")
  }
  fs <- fsApply(fs, function(ff){
    idx <- sample.int(nrow(ff), min(n, nrow(ff)))
    ff[idx,]
  })
  return(fs)
}
