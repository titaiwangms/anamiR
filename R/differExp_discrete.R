#' Find differential expression genes or miRNAs from given expression data
#'
#' This function will apply one of three statistical methods, including t.test,
#' wilcox.test and limma, to find differential expression genes or miRNAs with,
#' discrete phenotype data, and then filter the genes or miRNAs (rows) which
#' have bigger p-value than cutoff.
#'
#' @seealso \code{\link[stats]{t.test}} for Student's t-Test;
#'    \code{\link[stats]{wilcox.test}} for Wilcoxon Rank Sum and Signed Rank
#'    Tests.
#'
#' @return data expression data in matrix format, with sample name in columns and
#'    gene symbol or miRNA name in rows.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'    for input format.
#' @param class string. Choose one features from all rows of phenotype data.
#' @param method statistical method for finding differential genes or miRNAs,
#'    including "t.test", "wilcox.test", "limma". Default is "t.test".
#' @param limma.trend logical, only matter when limma is chosen to be the method.
#'    From function \code{\link[limma]{eBayes}}.
#' @param t_test.var logical, only matter when limma is chosen to be the method.
#'    Whether to treat the two variances as being equal. From function
#'    \code{\link[stats]{t.test}}
#' @param log2 logical, if this data hasn't been log2 transformed yet, this one
#'    should be TRUE. Default is FALSE.
#' @param p_value.cutoff an numeric value indicating a threshold of p-value
#'    for every genes or miRNAs (rows). Default is 0.05.
#' @param p_adjust.method Correction method for multiple testing. (If you are
#'    using DESeq for method, this param would not affect the result) From
#'    function \code{\link[stats]{p.adjust}}. Default is "BH".
#' @param foldchange an numeric value indicating a threshold of foldchange (log2)
#'    for every genes or miRNAs (rows). Default is 0.5.
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
#' ## Finding differential miRNA from miRNA expression data with t.test
#' mirna_d <- differExp_discrete(
#'    se = mirna_se,
#'    class = "ER",
#'    method = "t.test"
#' )
#'
#' @import stats
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment assays
#' @importFrom limma lmFit
#' @importFrom limma makeContrasts
#' @importFrom limma contrasts.fit
#' @importFrom limma eBayes
#' @importFrom DESeq2 DESeqDataSet
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#' @export
differExp_discrete <- function(
  se,
  class,
  method = c("t.test",
             "limma",
             "wilcox.test",
             "DESeq"),
  limma.trend = FALSE,
  t_test.var = FALSE,
  log2 = FALSE,
  p_value.cutoff = 0.05,
  p_adjust.method = "BH",
  foldchange = 0.5
) {

  data <- SummarizedExperiment::assays(se)[[1]]

  if (log2 %in% TRUE) {
    data <- log2(data)
  }

  method <- match.arg(method)
  pheno_data <- SummarizedExperiment::colData(se)[[class]]

  if (!is.null(pheno_data)) {
    pheno_data <- t(pheno_data)
    # seperate group
    if (method == "limma") {
      group <- t(pheno_data)
      type <- as.character(unique(unlist(group)))
      levels(group)[levels(group) == type[1]] <- 0
      levels(group)[levels(group) == type[2]] <- 1
    }
    group_1 <- which(pheno_data == levels(pheno_data)[1])
    group_2 <- which(pheno_data == levels(pheno_data)[2])

    #Fold Changes
    foldchange_cal <- function(da, gp1, gp2) {
      mean_gp1 <- apply(da[, gp1], 1, mean)
      mean_gp2 <- apply(da[, gp2], 1, mean)
      FC <- mean_gp1 - mean_gp2
      return(FC)
    }

    if (method != "DESeq") {
      FC <- foldchange_cal(data, group_1, group_2)
      p_value <- vector(mode = "numeric", length = nrow(data))
    }

    # t.test
    if (method == "t.test") {
      t_test <- function (da, gp1, gp2) {
        stats::t.test(da[gp1], da[gp2], var.equal = t_test.var)[["p.value"]]
      }
      p_value <- apply(data, 1, t_test, group_1, group_2)
    }

    # wilcoxon
    if (method == "wilcox.test") {
      wilcoxon <- function (da, gp1, gp2) {
        stats::wilcox.test(da[gp1], da[gp2])[["p.value"]]
      }
      p_value <- apply(data, 1, wilcoxon, group_1, group_2)
    }

    # limma  (trend)
    if (method == "limma") {
      design <- stats::model.matrix(~ 0 + group)
      fit <- limma::lmFit(data, design)
      mc <- limma::makeContrasts("group0 - group1", levels = design)
      fit2 <- limma::contrasts.fit(fit, mc)
      eb <- limma::eBayes(fit2, trend = limma.trend)
      p_value <- eb[["p.value"]]
    }

    #DESeq
    if (method == "DESeq") {
      tmp <- as.formula(paste("~", class))
      dds <- DESeq2::DESeqDataSet(se, design = tmp)
      dds <- DESeq2::DESeq(dds)
      res <- DESeq2::results(dds)
      p_value <- res[["pvalue"]]
      p_adjust <- res[["padj"]]
      FC <- res[["log2FoldChange"]]
    }

    # output
    if (method != "DESeq") {
      p_adjust <- stats::p.adjust(p = p_value, method = p_adjust.method)
    }
    idx <- which(p_adjust < p_value.cutoff)
    DE_data <- data[idx, ]
    DE_data <- cbind(DE_data, FC[idx], p_value[idx], p_adjust[idx])
    len_col <- ncol(DE_data)
    colnames(DE_data)[(len_col - 2):len_col] <- c("Fold-Change",
                                                   "P-Value",
                                                   "P-adjust")
    FC_rows <- abs(DE_data[, len_col - 2])
    DE_data <- DE_data[which(FC_rows > foldchange), ]
    return(DE_data)
  }
}
