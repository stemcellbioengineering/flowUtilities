

#' Estimate the logicle transform parameters for a GatingSet. Uses all the samples
#' within the GatingSet and returns a [flowWorkspace::transformerList] object that
#' can be applied to the GatingSet using [flowWorkspace::transform].
#' 
#' @param gs A GatingSet object
#' @param channels Character vector of channels for which the transformation is to be estimated
#' @param m The full width of the transformed display in asymptotic decades. Should be greater than zero
#' 
#' @returns A [flowWorkspace::transformerList] object
#' 
#' @export
estimateLogicleGS <- function(gs, channels, m = 4.5){
  if (!inherits(gs, "GatingSet")){
    stop("Argument 'gs' must be a GatingSet object")
  }
  # Merge into a single sample (GatingHierarchy)
  # Using estimateLogicle with flowFrame will raise error when added to GatingSet
  gh <- as(gs,"GatingHierarchy")
  
  # Estimate transform
  transList <- flowWorkspace::estimateLogicle(gh, channels, m = m)

  return (transList)
}