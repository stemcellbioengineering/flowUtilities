

#' Estimate the logicle transform parameters for a GatingSet or GatingHierarchy.
#' Uses all the samples within the GatingSet and returns a [flowWorkspace::transformerList] object that
#' can be applied using [flowWorkspace::transform].
#'
#' @param gs A GatingSet or GatingHierarchy object
#' @param channels Character vector of channels for which the transformation is to be estimated
#' @param m The full width of the transformed display in asymptotic decades. Should be greater than zero
#' @returns A [flowWorkspace::transformerList] object
#' @import flowWorkspace
#' @importFrom methods as
#' @export
estimateLogicleTransform <- function(gs, channels, m = 4.5){
  if (inherits(gs, "GatingSet")){
    # Merge into a single sample (GatingHierarchy)
    # Using estimateLogicle with flowFrame will raise error when added to GatingSet
    gh <- as(gs,"GatingHierarchy")
  }else if (inherits(gs, "GatingHierarchy")){
    gh <- gs
  }else{
    stop("Argument 'gs' must be a GatingSet or GatingHierarchy object")
  }
  # Estimate transform
  transList <- estimateLogicle(gh, channels, m = m)

  return (transList)
}
