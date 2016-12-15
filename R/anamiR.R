#' anamiR: An integrated analysis package of miRNA and mRNA expression
#' data.
#'
#' The anamiR package is used to identify miRNA-target genes interactions.
#' The anamiR package provides a whole workflow, which contains important
#' functions: `normalization`, `differExp_discrete`, `negative_cor`,
#' `miR_converter`, `database_support`, `enrichment`.
#'
#' @section normalization: The normalization function is used to
#'   normalize the expression data with one of three methods,
#'   including normal, quantile, rank.invariant.
#'
#' @section differExp_discrete: The differExp_discrete function
#'   is used to find the differential genes or miRNAs from given
#'   expression data with one of three statistical methods,
#'   including t.test, wilcox.test,limma and DESeq. The miRNA
#'   would remain if its p-value lower than the cutoff value.

#' @section miR_converter: The miR_annotation function is used
#'   to convert the older miRNA annotation to the miRBase 21
#'   version.
#'
#' @section negative_cor: The negative_cor  function is used to
#'   identify the possible miRNA-target gene interactions from
#'   given miRNA and mRNA expression data by caculating the
#'   correlation coefficient between each miRNA and gene.
#'   interaction would remain if its correlation coefficient
#'   is negative and lower than cutoff value.
#'
#' @section database_support: The database_support function would
#'   search information about miRNA-target gene interactions from
#'   an integrated database, which contains 8 algorithm predicted
#'   databases and 2 experiment validated databases. Eventually
#'   return a big table, which is in data.frame format and contains
#'   extra 10 columns for those 10 databases to count if interactions
#'   were predicted or validated by these databases.
#'
#' @section enrichment: The enrichment function is used to do the
#'   functional analysis from the output of `database_support`.
#'   Not only p-value from hypergeometric test but empirical p-value
#'   from 10000 times of permutation would be provided by this function.
#'
#' @docType package
#' @name anamiR
NULL
