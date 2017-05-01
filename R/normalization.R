#' Normalize expression data
#'
#' This function will normalize the given expression data and return
#' it in the same data format.
#'
#' @seealso \code{\link[limma]{normalizeQuantiles}} for quantile normalization;
#'    \code{\link[lumi]{rankinvariant}} for rank invariant normalization.
#'
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'    for return object.
#'
#' @param data expression data in matrix format, with sample name in columns
#'    and gene symbol or miRNA name in rows.
#' @param method normalization methods, including "quantile", "normal",
#'    "rank.invariant". Default is "quantile". As for method "normal",
#'    we trim the extreme value and calculate the mean in the data.
#'
#' @examples
#' ## Use the internal dataset
#' data("mirna", package = "anamiR", envir = environment())
#'
#' ## Normalize miRNA expression data
#' normalization(data = mirna, method = "quantile")
#'
#' @importFrom lumi rankinvariant
#' @import limma
#' @export
normalization <- function(
  data,
  method = c("quantile",
             "normal",
             "rank.invariant")
) {

  if (method %in% "normal") {
    data_trm <- apply(data, 2, mean, trim = 0.02)
    sd <- apply(data, 2, sd)
    median <- apply(data, 2, median)
    data_trm_mean <- mean(data_trm)
    data_trmean <- data / data_trm * data_trm_mean
    data <- data_trmean
  }

  if (method %in% "quantile") {
    data <- limma::normalizeQuantiles(data)
  }

  if (method %in% "rank.invariant") {
    data <- lumi::rankinvariant(as.matrix(data))
    data <- as.data.frame(data)
  }

  return(data)
}
