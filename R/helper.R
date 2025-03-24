
#' evaluate signal detection error
#'
#' @param true_signal_indices vector, contains the indices of true signals
#' @param estimated_signal_indices vector, contains the indices of estimated signals
#'
#' @returns
#'  \describe{
#'    \item{type_I_error}{Type-I error of signal detection}
#'    \item{type_II_error}{Type-II error of signal detection}
#'  }
#' @export
signal_detection_error <- function(true_signal_indices, estimated_signal_indices) {
  true_positive <- length(intersect(true_signal_indices, estimated_signal_indices))
  false_positive <- length(setdiff(estimated_signal_indices, true_signal_indices))
  false_negative <- length(setdiff(true_signal_indices, estimated_signal_indices))

  type_I_error <- false_positive / length(estimated_signal_indices)
  type_II_error <- false_negative / length(true_signal_indices)

  return(list(type_I_error = type_I_error, type_II_error = type_II_error))
}

#' evaluate estimation error
#'
#' @param true_pc true principal component
#' @param estimated_pc estimated principal component
#'
#' @returns
#'  \describe{
#'    \item{alignment}{a float between 0 and 1, alignment between true and estimated PC}
#'  }
#' @export
estimation_error <- function(true_pc, estimated_pc){

  alignment <- abs(sum(true_pc * estimated_pc))

  return(alignment)
}
