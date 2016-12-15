#' Find differential expression genes or miRNAs from given expression data
#'
#' This function will apply linear regression model to find differential expression
#' genes or miRNAs with continuous phenotype data,and then filter the genes or
#' miRNAs (rows) which have bigger p-value than cutoff.
#'
#' @seealso \code{\link[stats]{lm}} for fitting linear models.
#'
#' @return data expression data in matrix format, with sample name in columns and
#'    gene symbol or miRNA name in rows.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'    for input format.
#' @param class string. Choose one features from all rows of phenotype data.
#' @param log2 logical, if this data hasn't been log2 transformed yet, this one
#'    should be TRUE. Default is FALSE.
#' @param p_value.cutoff an numeric value indicating a threshold of p-value
#'    for every genes or miRNAs (rows). Default is 0.05.
#'
#' @examples
#' ## Use the internal dataset
#' data("mirna", package = "anamiR", envir = environment())
#' data("pheno.mirna", package = "anamiR", envir = environment())
#'
#' ## SummarizedExperiment class
#' require(SummarizedExperiment)
#' mirna_se <- SummarizedExperiment(
#'  assays = SimpleList(counts=mirna),
#'  colData = pheno.mirna)
#'
#' ## Finding differential miRNA from miRNA expression data with lm
#' differExp_continuous(
#'     se = mirna_se, class = "Survival"
#' )
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment assays
#' @import stats
#' @export
differExp_continuous <- function(
  se,
  class,
  log2 = FALSE,
  p_value.cutoff = 0.05
) {

  if (log2 %in% TRUE) {
    data <- log2(data)
  }

  data <- SummarizedExperiment::assays(se)[[1]]
  pheno_data <- SummarizedExperiment::colData(se)[[class]]

  if (!is.null(pheno_data)) {
  # linear model
    p_value <- vector(mode = "numeric", length = nrow(data))
    a <- stats::lm(t(data) ~ pheno_data)
    a_sum <- summary(a)
    for (i in seq_len(length(a_sum))) {
      p_value[i] <- a_sum[[i]][["coefficients"]][2, 4]
    }
  }
  # output
  idx <- which(p_value < p_value.cutoff)
  DE_data <- data[idx, ]
  DE_data <- cbind(DE_data, p_value[idx])
  len_col <- ncol(DE_data)
  colnames(DE_data)[len_col] <- c("P-Value")
  return(DE_data)
}
