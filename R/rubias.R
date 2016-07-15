#' rubias: A package for bias correction in hierarchical GSI
#'
#' Blah blah blah.
#'
#' @section rubias functions:
#' The rubias functions ...
#'
#' @section example data:
#' \code{alewife} & \code{blueback} are examples of good times to use rubias; there are
#' clearly observable biases based on the number of collections in a reporting unit.
#'  \code{chinook} is an example of a dataset in which examination of the graphs from
#' \code{bias_comparison} reveal that rubias is not appropriate; there was no appreciable
#' bias to begin with, and so both the misassignment-scaling and parametric boostrap methods
#' have no consistent effect on the mean residual, and often increase mean squared error
#'
#' @docType package
#' @name rubias
#' @importFrom Rcpp evalCpp
#' @useDynLib rubias
NULL




