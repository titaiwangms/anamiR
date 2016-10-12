#' Using Fold-Chang information to draw a heatmap
#'
#' This function would base on Fold-Change information
#' from the output of \link{database_support}, and
#' show a plot to users. Note that if miRNA-gene
#' interactions (row) from input are larger than 100,
#' the lable in plot would be unclear.
#'
#' @return heatmap plot with Fold-change of miRNA and
#'    gene in column names and name of miRNA and gene
#'    in row names.
#'
#' @seealso \code{\link[stats]{heatmap}} for plot.
#'
#' @param sup_data matrix format generated from
#'    \link{database_support}.
#'
#' @examples
#' ## Use the internal dataset
#' data("mirna", package = "anamiR", envir = environment())
#' data("pheno.mirna", package = "anamiR", envir = environment())
#' data("mrna", package = "anamiR", envir = environment())
#' data("pheno.mrna", package = "anamiR", envir = environment())
#'
#' ## SummarizedExperiment class
#' require(SummarizedExperiment)
#' mirna_se <- SummarizedExperiment(
#'  assays = SimpleList(counts=mirna),
#'  colData = pheno.mirna)
#'
#' ## SummarizedExperiment class
#' require(SummarizedExperiment)
#' mrna_se <- SummarizedExperiment(
#'  assays = SimpleList(counts=mrna),
#'  colData = pheno.mrna)
#'
#' ## Finding differential miRNA from miRNA expression data with t.test
#' mirna_d <- differExp_discrete(
#'    se = mirna_se,
#'    class = "ER",
#'    method = "t.test"
#' )
#'
#' ## Finding differential mRNA from mRNA expression data with t.test
#' mrna_d <- differExp_discrete(
#'    se = mrna_se,
#'    class = "ER",
#'    method = "t.test"
#' )
#'
#' ## Convert annotation to miRBse 21
#' mirna_21 <- miR_converter(data = mirna_d, original_version = 17)
#'
#' ## Correlation
#' cor <- negative_cor(mrna_data = mrna_d, mirna_data = mirna_21)
#'
#' ## Intersect with known databases
#' sup <- database_support(cor_data = cor)
#'
#' ## Draw heatmap
#' heat_vis(sup)
#'
#' @import stats
#' @export

heat_vis <- function(
  sup_data
) {
  fc_mirna <- as.numeric(sup_data[, 18])
  fc_gene <- as.numeric(sup_data[, 20])
  fc_mat <- cbind(fc_mirna, fc_gene)
  row.names(fc_mat) <- paste(sup_data[, 1],"-",sup_data[, 2])
  colnames(fc_mat) <- c(colnames(sup_data)[18],
                        colnames(sup_data)[20])
  stats::heatmap(fc_mat, cexRow = 0.8, cexCol = 0.8)
}
