#' The 'bbmix' package.
#'
#' @description Bayesian Beta-Binomial mixture model for RNA-seq genotyping
#'
#' @docType package
#' @name bbmix-package
#' @aliases bbmix
#' @useDynLib bbmix, .registration = TRUE
#' @import R.utils
#' @import data.table
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom stats complete.cases rbeta rgamma sd
#' @importFrom utils write.table
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.18.2. https://mc-stan.org
#'
NULL
