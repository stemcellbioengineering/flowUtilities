

#' Constructs a singlet gate by applying a robust linear model.
#'
#' A convenience function for [flowStats::singletGate] that constructs a gate for
#' each sample within a GatingSet, flowSet, or cytoset and returns a list of gates.
#' Additionally, when a GatingSet is provided, the singlet gate can be constructed after filtering by an existing gate.
#'
#' @param gs A GatingSet, flowSet, or cytoset object
#' @param area Character giving the channel name that records the signal intensity as peak area
#' @param height Character giving the channel name that records the signal intensity as peak height
#' @param sidescatter UNUSED Character giving an optional channel name for the sidescatter signal (by default, ignored)
#' @param prediction_level A numeric value between 0 and 1 specifying the level to use for the prediction bands
#' @param subsample_pct A numeric value between 0 and 1 indicating the percentage of observations that should be randomly selected from x to construct the gate (by default, no subsampling is performed)
#' @param wider_gate Logical value. If TRUE, the prediction bands used to construct the singlet gate use the robust fitted weights, which increase prediction uncertainty, especially for large FSC-A. This leads to wider gates, which are sometimes desired
#' @param parentId Character specify the parent node name used to filter the data. Only used when a GatingSet is provided (by default, the 'root' node, no filtering is performed)
#' @param gateId Character specifying the name for the gate that is returned
#' @param maxit the limit on the number of IWLS iterations passed to [MASS::rlm]
#' @param ... Additional arguments passed to [MASS::rlm]
#' @returns Returns a list of [flowCore::polygonGate] objects to use for filtering
#' @examples
#' \dontrun{
#' # Construct gates after filtering by existing "NonDebris" gate
#' gates <- autoSingletGate(gs,
#'                          area = "FSC-A",
#'                          height = "FSC-H",
#'                          parentId = "NonDebris",
#'                          gateId = "Singlets")
#'
#' # Add gate to GatingSet
#' gs_pop_add(gs, gates, parent = "NonDebris")
#' }
#' @importFrom flowStats singletGate
#' @import flowCore
#' @import flowWorkspace
#' @export
autoSingletGate <- function(gs,
                          area="FSC-A",
                          height="FSC-H",
                          sidescatter=NULL,
                          prediction_level=0.99,
                          subsample_pct=NULL,
                          wider_gate=FALSE,
                          parentId="root",
                          gateId="singlets",
                          maxit=10,
                          ...){
  # Validate inputs
  if (inherits(gs, c("flowframe","cytoframe"))) {
    stop("This function is for constructing gates for multiple samples. For flowFrames or cytoframes, use flowStats::singletGate directly")
  }else if (!inherits(gs, c("GatingSet","GatingHierarchy","flowSet","cytoset"))) {
    stop("Argument 'gs' must be a GatingSet, GatingHierarchy, flowSet or cytoset object")
  }else if (inherits(gs, c("GatingSet","GatingHierarchy"))){
    # Get data from GatingSet after filtering by parentId
    # Returns a cytoset
    gs <- gs_pop_get_data(gs, parentId)
  }

  if (!(area %in% colnames(gs)) | !(height %in% colnames(gs))){
    stop(sprintf("One of channels %s or %s not in dataset", area, height))
  }

  # Estimate per sample using RLM
  gates <- fsApply(gs, function(cf) {singletGate(cf,
                                     area=area,
                                     height=height,
                                     #sidescatter=sidescatter,
                                     prediction_level=prediction_level,
                                     subsample_pct=subsample_pct,
                                     wider_gate=wider_gate,
                                     filterId=gateId,
                                     maxit=maxit,
                                     ...
                                     )})
  return (gates)
}

#' Automated gating of elliptical cell populations in 2D.
#'
#' A convenience function for [flowStats::lymphGate] that constructs a gate for
#' each sample within a GatingSet, flowSet, or cytoset and returns a list of gates.
#' Additionally, when a GatingSet is provided, the gate can be constructed after
#' filtering by an existing gate.
#'
#' @param gs A GatingSet, flowSet, or cytoset object
#' @param x Character providing a channel name
#' @param y Character providing a channel name
#' @param preselection A character giving a channel name in gs, or a named list of numerics specifying the initial rough preselection. The latter gets passed on to [flowCore::rectangleGate]. If not provided, fits a regular [flowStats::norm2Filter] (the default)
#' @param scale The scaleFactor parameter that gets passed on to [flowStats::norm2Filter]
#' @param bwFac The bandwidth factor that gets passed on to [flowStats::curv1Filter]
#' @param parentId Character specify the parent node name used to filter the data. Only used when a GatingSet is provided (by default, the 'root' node, no filtering is performed)
#' @param gateId Character specifying the name for the gate that is returned
#' @param ... Additional arguments passed to [flowStats::lymphGate]
#' @returns Returns a list of [flowCore::ellipsoidGate] objects to use for filtering
#' @examples
#' \dontrun{
#' # Construct gates after filtering by existing "Singlets" gate
#' gates <- autoLymphGate(gs,
#'                        x = "FSC-A",
#'                        y = "SSC-A",
#'                        parentId = "Singlets",
#'                        gateId = "Lymphocytes")
#'
#' # Add gate to GatingSet
#' gs_pop_add(gs, gates, parent = "Singlets")
#' }
#' @importFrom flowStats lymphGate
#' @import flowCore
#' @import flowWorkspace
#' @export
autoLymphGate <- function(gs,
                          x="FSC-A",
                          y="SSC-A",
                          preselection=NULL,
                          scale=2.5,
                          bwFac=1.3,
                          parentId="root",
                          gateId="lymphocytes",
                          ...){
  # Validate inputs
  if (inherits(gs, c("flowframe","cytoframe"))) {
    stop("This function is for constructing gates for multiple samples. For flowFrames or cytoframes, use flowStats::lymphtGate directly")
  }else if (!inherits(gs, c("GatingSet","GatingHierarchy","flowSet","cytoset"))) {
    stop("Argument 'gs' must be a GatingSet, GatingHierarchy, flowSet or cytoset object")
  }else if (inherits(gs, c("GatingSet","GatingHierarchy"))){
    # Get data from GatingSet after filtering by parentId
    # Returns a cytoset
    gs <- gs_pop_get_data(gs, parentId)
  }

  if (!(x %in% colnames(gs)) | !(y %in% colnames(gs))){
    stop(sprintf("One of channels %s or %s not in dataset", x, y))
  }

  # Estimate per sample using lymphGate
  gates <- lymphGate(gs,
                     channels = c(x,y),
                     preselection = preselection,
                     scale = scale,
                     bwFac = bwFac,
                     filterId=gateId,
                     ...)

  return (gates)
}

#' Constructs a quadrant gate after estimating the quantiles from controls.
#' Requires the phenoData to contain a column 'ctrlId', containing names of control
#' samples that match the channel names. This allows the function to select the correct
#' control for calculating the quantile. The controls are usually fluorescent minus one (FMO).
#'
#' @param gs A GatingSet object
#' @param x Character providing a channel name
#' @param y Character providing a channel name
#' @param probs Numeric probability used to estimate the quantile (default: 0.99)
#' @param ctrlId Character name of column in phenoData matching control samples to channel names (default: "fmo")
#' @param parentId Character specify the parent node name used to filter the data (by default, the 'root' node, no filtering is performed)
#' @param gateId Character specifying the name for the gate that is returned
#' @returns A [flowCore::quadGate] object
#' @examples
#' \dontrun{
#' # Construct gates after filtering by existing "Live" gate
#' gate <- autoQuadGate(gs, x = "CD4", y = "CD8", parentId = "Live")
#'
#' # A quadGate becomes four rectangleGate when added to a GatingSet
#' # The names will be created automatically or, we can specify the names
#' # when adding to GatingSet (order is clockwise from upper left)
#' names <- c("CD4-CD8+","CD4+CD8+","CD4+CD8-","CD4-CD8-")
#'
#' # Add gate to GatingSet
#' gs_pop_add(gs, gate, parent = "Live", names = names)
#' }
#' @import flowCore
#' @import flowWorkspace
#' @export
autoQuadGate <- function(gs,
                        x,
                        y,
                        probs = 0.99,
                        ctrlId = "fmo",
                        parentId="root",
                        gateId="quadGate"){
  # Validate inputs
  if (inherits(gs, c("flowframe","cytoframe"))) {
    stop("This function is for constructing gates for multiple samples")
  }else if (!inherits(gs, c("GatingSet","GatingHierarchy","flowSet","cytoset"))) {
    stop("Argument 'gs' must be a GatingSet, GatingHierarchy, flowSet or cytoset object")
  }else if (inherits(gs, c("GatingSet","GatingHierarchy"))){
    # Get data from GatingSet after filtering by parentId
    # Returns a cytoset
    gs <- gs_pop_get_data(gs, parentId)
  }
  # Check for matching names
  if (!(x %in% colnames(gs)) | !(y %in% colnames(gs))){
    stop(sprintf("Argument 'x': %s or 'y': %s not in dataset", x, y))
  }else if (!(ctrlId %in% colnames(pData(gs)))){
    stop(sprintf("Argument 'ctrlId': %s must be a column in phenoData", ctrlId))
  }else if (!(x %in% pData(gs)[[ctrlId]]) | !(y %in% pData(gs)[[ctrlId]])){
    stop(sprintf("Argument 'x': %s or 'y': %s not in phenoData column %s", x, y, ctrlId))
  }

  # Get quantile for x
  gs_x <- gs[pData(gs)[[ctrlId]]==x]
  qt_x <- calculate_quantile(gs_x,
                             channels = x,
                             probs = probs)
  # Get quantile for y
  gs_y <- gs[pData(gs)[[ctrlId]]==y]
  qt_y <- calculate_quantile(gs_y,
                             channels = y,
                             probs = probs)
  # Merge quantiles
  quantiles <- c(qt_x, qt_y)
  print(quantiles)
  # Build gate
  gate <- quadGate(quantiles, filterId = gateId)

  return (gate)
}

#' Constructs a 6 vertex diagonal polygon gate. Useful for gating singlets using
#' FSC-A x FSC-H or gating CD3+TCRab+ T-cells.
#'
#' The polygon has a maximum width defined by the width argument, a maximum height
#' defined by the height argument, and a width between diagonal lines defined by diag_width.
#' The bottom left corner is positioned at the origin argument.
#'
#' @param x The name of the channel along the width (x-axis)
#' @param y The name of the channel along the height (y-axis)
#' @param width The total width of the polygon
#' @param height The total height of the polygon
#' @param origin A numeric vector of length 2 defining the position of the bottom left corner of the gate (the default is c(0,0))
#' @param diag_width The width between the two diagonal lines (the default is 0.5)
#' @param gateId Character specifying the name for the gate that is returned
#' @returns Returns a [flowCore::polygonGate] object to use for filtering
#' @examples
#' \dontrun{
#' # Construct a singlets gate
#' gates <- diagPolygonGate(gs,
#'                          x = "FSC-A",
#'                          y = "FSC-H",
#'                          width = 3e6,
#'                          height = 3.5e6,
#'                          diag_width = 7e5,
#'                          gateId = "Singlets")
#'
#' # Add gate to GatingSet
#' gs_pop_add(gs, gates, parent = "root")
#' }
#' @import flowCore
#' @export
diagPolygonGate <- function(x,
                            y,
                            width,
                            height,
                            origin = c(0,0),
                            diag_width = 0.5,
                            gateId = "PolygonGate"){
  # Input validation
  if (!inherits(origin,"numeric")){
    stop("Origin must be a numeric vector")
  }else if(!length(origin) == 2){
    stop("Origin must have length 2")
  }

  # Calculate angle of diagonal
  theta <- atan(height/width)
  # Length of each base
  base_height <- diag_width * sin(pi/2 - theta)
  base_width <- diag_width * cos(pi/2 - theta)
  # Make vertices for x,y
  x_vert <- c(0, 0, width - base_width, width, width, base_width)
  y_vert <- c(0,base_height, height, height, height - base_width, 0)
  # Add origin offset
  x_vert = x_vert + origin[[1]]
  y_vert = y_vert + origin[[2]]
  # Combine and name
  vertices <- list(x_vert,y_vert)
  names(vertices) <- c(x, y)
  # Construct gate
  gate <- polygonGate(vertices, filterId = gateId)
  return(gate)
}
