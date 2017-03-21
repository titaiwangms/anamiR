#' This function will do GSEA analysis.
#'
#' This function will do GSEA analysis through the function
#' \code{\link[gage]{gage}}. After obtaining the ranking of
#' pathways, this function will choose the top five (default)
#' pathaways, and then find the related miRNAs based on their
#' gene set.
#'
#' @return list format containing both selected gene and miRNA
#'    expression data for each chosen pathway.
#'
#' @seealso \code{\link[gage]{gage}} for GSEA analysis.
#'
#' @param mrna_se \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'    for input format and it contains mRNA information.
#' @param mirna_se \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'    for input format, and it contains miRNA information.
#' @param class string. Choose one features from all rows of phenotype data.
#' @param eg2sym logical. conversion between Entrez Gene IDs
#'    and official gene symbols for human genes.
#' @param compare character, if the length of case is the same
#'    as control, use "paired".Default is "unpaired".
#' @param pathway_num The number of chosen pathways from the
#'    result of GSEA analysis.
#'
#' @examples
#'
#' require(data.table)
#'
#' ## Load example data
#' aa <- system.file("extdata", "GSE19536_mrna.csv", package = "anamiR")
#' mrna <- fread(aa, fill = TRUE, header = TRUE)
#'
#' bb <- system.file("extdata", "GSE19536_mirna.csv", package = "anamiR")
#' mirna <- fread(bb, fill = TRUE, header = TRUE)
#'
#' cc <- system.file("extdata", "pheno_data.csv", package = "anamiR")
#' pheno.data <- fread(cc, fill = TRUE, header = TRUE)
#'
#' ## adjust data format
#' mirna_name <- mirna[["miRNA"]]
#' mrna_name <- mrna[["Gene"]]
#' mirna <- mirna[, -1]
#' mrna <- mrna[, -1]
#' mirna <- data.matrix(mirna)
#' mrna <- data.matrix(mrna)
#' row.names(mirna) <- mirna_name
#' row.names(mrna) <- mrna_name
#' pheno_name <- pheno.data[["Sample"]]
#' pheno.data <- pheno.data[, -1]
#' pheno.data <- as.matrix(pheno.data)
#' row.names(pheno.data) <- pheno_name
#'
#' ## SummarizedExperiment class
#' require(SummarizedExperiment)
#' mirna_se <- SummarizedExperiment(
#'  assays = SimpleList(counts=mirna),
#'  colData = pheno.data)
#'
#' mrna_se <- SummarizedExperiment(
#'  assays = SimpleList(counts=mrna),
#'  colData = pheno.data)
#'
#' table <- GSEA_ana(mrna_se = mrna_se,
#'  mirna_se = mirna_se, class = "ER",
#'  pathway_num = 2)
#'
#' @import gage
#' @import RMySQL
#' @import DBI
#' @importFrom S4Vectors head
#' @export
GSEA_ana <- function(
  mrna_se,
  mirna_se,
  class,
  compare = "unpaired",
  eg2sym = TRUE,
  pathway_num = 5
) {

  mrna <- SummarizedExperiment::assays(mrna_se)[[1]]
  mirna <- SummarizedExperiment::assays(mirna_se)[[1]]
  pheno_data <- SummarizedExperiment::colData(mrna_se)[[class]]

  if (eg2sym %in% TRUE ) {
    row.names(mrna) <- gage::sym2eg(row.names(mrna))
  }

  na_row <- which(is.na(row.names(mrna)))
  mrna <- mrna[-na_row, ]

  mrna_name <- gage::eg2sym(row.names(mrna))
  mirna_name <- row.names(mirna)

  pheno_data <- t(pheno_data)
  # seperate group
  gp1 <- which(pheno_data == levels(pheno_data)[1])
  gp2 <- which(pheno_data == levels(pheno_data)[2])

  kegg.p <- gage::gage(mrna, gsets = kegg.gs,
                       ref = gp2, samp = gp1, compare = compare)

  pathway <- S4Vectors::head(kegg.p$greater[, 1:pathway_num],
                  pathway_num)

  # connect with db
  db <- RMySQL::dbConnect(RMySQL::MySQL(), user = "visitor",
                          password = "visitor",
                          dbname = "visitor",
                          host = "anamir.cgm.ntu.edu.tw")
  gene_list <- c()
  mirna_list <- c()
  table <- list()

  for (i in seq_len(pathway_num)) {
    if (grepl("\\(", row.names(pathway)[i])) {
      pat <- strsplit(x = row.names(pathway)[i],
                      split = "\\(")[[1]][1]
      number <- grep(pat, names(kegg.gs))
    } else {
      number <- grep(row.names(pathway)[i], names(kegg.gs))
    }
    gset_id <- kegg.gs[[number]]
    gset_sym <- gage::eg2sym(gset_id)
    for (j in seq_len(length(gset_sym))) {
      gene <- gset_sym[j]
      query <- paste0("SELECT `miRNA_21`, `gene_symbol`
                      FROM `all_hsa` where gene_symbol like '",
                      gene, "' ;")
      tmp <- DBI::dbGetQuery(db, query)
      row_need <- which(tmp[["miRNA_21"]] %in% mirna_name ||
                          tmp[["gene_symbol"]] %in% mrna_name)
      tmp <- tmp[row_need, ]
      gene_list <- c(gene_list, tmp[["gene_symbol"]])
      mirna_list <- c(mirna_list, tmp[["miRNA_21"]])
    }
    mirna_list <- unique(mirna_list)
    gene_list <- unique(gene_list)

    mirow <- which(mirna_name %in% mirna_list)
    mrow <- which(mrna_name %in% gene_list)

    mirna_tmp <- mirna[mirow, ]
    mrna_tmp <- mrna[mrow, ]
    row.names(mrna_tmp) <- gage::eg2sym(row.names(mrna_tmp))

    table[[2 * i - 1]] <- mirna_tmp
    names(table)[[2 * i - 1]] <- paste(row.names(pathway)[i],
                                       "- mirna")
    table[[2 * i]] <- mrna_tmp
    names(table)[[2 * i]] <- paste(row.names(pathway)[i],
                                       "- mrna")
  }

  #disconnect db
  cons <- RMySQL::dbListConnections(RMySQL::MySQL())
  for (con in cons) RMySQL::dbDisconnect(con)

  return(table)
}
