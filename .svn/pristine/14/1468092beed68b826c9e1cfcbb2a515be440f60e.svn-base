#' Pipeline of anamiR is applied to given output from \link{GSEA_ana}.
#'
#' This function will use \link{differExp_discrete} and
#' \link{negative_cor} to do the deeper analysis of given data which
#' is from \link{GSEA_ana}.
#'
#' @return list format containing matrix for each chosen pathway.
#'    The format of matrix is like the output from \link{negative_cor}.
#'
#' @seealso \link{differExp_discrete} and \link{negative_cor}.
#'
#' @param table list format containing both selected gene and miRNA
#'    expression data for each chosen pathway. output of \link{GSEA_ana}
#' @param pheno.data phenotype data.
#' @param class string. Choose one features from all rows of phenotype data.
#' @param DE_method statistical method for finding differential genes or miRNAs,
#'    including "t.test", "wilcox.test", "limma". Default is "t.test".
#' @param limma.trend logical, only matter when limma is chosen to be the method.
#'    From function \code{\link[limma]{eBayes}}.
#' @param t_test.var logical, only matter when limma is chosen to be the method.
#'    Whether to treat the two variances as being equal. From function
#'    \code{\link[stats]{t.test}}
#' @param log2 logical, if this data hasn't been log2 transformed yet, this one
#'    should be TRUE. Default is FALSE.
#' @param p_adjust.method Correction method for multiple testing. (If you are
#'    using DESeq for method, this param would not affect the result) From
#'    function \code{\link[stats]{p.adjust}}. Default is "BH".
#' @param cor_cut an numeric value indicating a threshold of correlation
#'    coefficient for every potential miRNA-genes interactions. Default is -0.3,
#'    however, if no interaction pass the threshold, this function would add
#'    0.2 value in threshold until at least one interaction passed the threshold.
#'
#' @examples
#' ## Load example data
#'
#' require(data.table)
#'
#' cc <- system.file("extdata", "pheno_data.csv", package = "anamiR")
#' pheno.data <- fread(cc, fill = TRUE, header = TRUE)
#'
#' ## adjust data format
#' pheno_name <- pheno.data[["Sample"]]
#' pheno.data <- pheno.data[, -1]
#' pheno.data <- as.matrix(pheno.data)
#' row.names(pheno.data) <- pheno_name
#' data(table)
#'
#' result <- GSEA_res(table = table, pheno.data = pheno.data,
#'  class = "ER", DE_method = "limma")
#'
#' @importFrom  SummarizedExperiment SummarizedExperiment
#' @importFrom  S4Vectors SimpleList
#' @export
GSEA_res <- function(
  table,
  pheno.data,
  class,
  DE_method = c("t.test",
             "limma",
             "wilcox.test",
             "DESeq"),
  limma.trend = FALSE,
  t_test.var = FALSE,
  log2 = FALSE,
  p_adjust.method = "BH",
  cor_cut = -0.3
) {
  cor_table <- list()
  #workflow
  n <- length(table) / 2
  for (i in seq_len(n)) {
    mirna_data <- table[[2 * i - 1]]
    mrna_data <- table[[2 * i]]

    mirna_se <- SummarizedExperiment::SummarizedExperiment(
      assays = S4Vectors::SimpleList(counts=mirna_data),
      colData = pheno.data)
    mrna_se <- SummarizedExperiment::SummarizedExperiment(
      assays = S4Vectors::SimpleList(counts=mrna_data),
      colData = pheno.data)

    mirna_dif <- differExp_discrete(se = mirna_se,
                                    class = class, method = DE_method,
                                    limma.trend = limma.trend,
                                    t_test.var = t_test.var,
                                    log2 = log2, logratio = 0,
                                    p_value.cutoff = 1,
                                    p_adjust.method = p_adjust.method)
    mrna_dif <- differExp_discrete(se = mrna_se,
                                    class = class, method = DE_method,
                                    limma.trend = limma.trend,
                                    t_test.var = t_test.var,
                                    log2 = log2, logratio = 0,
                                    p_value.cutoff = 1,
                                    p_adjust.method = p_adjust.method)

    cor <- negative_cor(mrna_data = mrna_dif, mirna_data = mirna_dif,
                        cut.off = cor_cut)
    cor_table[[i]] <- cor
    name <- strsplit(names(table)[2 * i], "-")[[1]][1]
    names(cor_table)[i] <- name
  }

  return(cor_table)
}
