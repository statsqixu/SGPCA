% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{signal_detection_error}
\alias{signal_detection_error}
\title{evaluate signal detection error}
\usage{
signal_detection_error(true_signal_indices, estimated_signal_indices)
}
\arguments{
\item{true_signal_indices}{vector, contains the indices of true signals}

\item{estimated_signal_indices}{vector, contains the indices of estimated signals}
}
\value{
\describe{
\item{type_I_error}{Type-I error of signal detection}
\item{type_II_error}{Type-II error of signal detection}
}
}
\description{
evaluate signal detection error
}
