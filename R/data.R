#' @title CRC subset data
#'
#' @description A small example dataset from the colorectal cancer (CRC) plasma proteome study.
#' It contains a pre-filtered expression matrix and sample group labels
#' (HC = healthy control, Pre = pre-treatment, Post = post-treatment).
#'
#' @usage data(CRC_subset)
#'
#' @format A list with two elements:
#' \describe{
#'   \item{expr}{A numeric matrix of features (rows) × samples (columns).}
#'   \item{group}{A character vector of sample groups, length = ncol(expr).}
#' }
#'
#' @details
#' The data are a subset of the discovery cohort used in the SwitchClass analysis.
#' Useful for testing or demonstrating functions such as \code{label_switch_classify()}.
#'
#' @source Internal preprocessed subset from \href{https://www.nature.com/articles/s41467-024-44911}{Li et al. 2024, Nat Commun}.
#'
#' @examples
#' data(CRC_subset)
#' dim(CRC_subset$expr)
#' table(CRC_subset$group)
#'
"CRC_subset"
