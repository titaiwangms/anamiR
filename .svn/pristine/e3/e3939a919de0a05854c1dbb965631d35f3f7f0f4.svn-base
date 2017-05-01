#' Convert miRNA annotation to the miRBase 21 version
#'
#' This function will convert the miRNA names from the data frame, which
#' is produced by \link{differExp_discrete}, to the miRBase 21 version
#' of miRNA annotation. If the input contains hundreds of miRNAs, it
#' would take a few minutes to convert all of them.
#'
#' @return expression data in data.frame format, with sample name in columns and
#'    miRNA name for miRBase version 21 in rows.
#'
#' @param data expression data in data.frame format, with sample name in columns and
#'    miRNA name in rows.
#' @param remove_old logical value, if the miRNA is deleted in miRBase 21, should
#'    it be removed from row? Default is TRUE.
#' @param original_version the original version of miRNA  in input matrix. This
#'    one is necessary.
#' @param latest_version choose an interger under 21, and this function would
#'    convert miRNA annotation to that version. Default is 21.
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
#' ## Convert annotation to miRBse 21
#' mirna_21 <- miR_converter(data = mirna_d, original_version = 17)
#'
#' @import RMySQL
#' @import DBI
#' @export
miR_converter <- function (
  data,
  remove_old = TRUE,
  original_version,
  latest_version = 21
) {
  mirna_list <- row.names(data)
  db <- RMySQL::dbConnect(RMySQL::MySQL(),
                          user = "visitor",
                          password = "visitor",
                          dbname = "visitor",
                          host = "anamir.cgm.ntu.edu.tw")
  mirna_21 <- c()
  for (mirna in mirna_list) {
    query <- paste0("SELECT * FROM `all_ver_miRNA` WHERE name like '",
                    mirna,
                    "';")
    annotation <- DBI::dbGetQuery(db, query)
    if (length(annotation[["version"]]) == 0) {
      mirna_21 <- c(mirna_21, "None")
    } else if (latest_version %in% annotation[["version"]]) {
      mirna_21 <- c(mirna_21, mirna)
    } else {
      accession <- annotation[["accession"]][1]
      query2 <- paste0("SELECT * FROM `all_ver_miRNA` WHERE accession like '",
                        accession,
                        "' and version like '",
                        latest_version,
                        "';")
      annotation2 <- DBI::dbGetQuery(db, query2)
      if (length(annotation2[["version"]]) == 0) {
        mirna_21 <- c(mirna_21, "None")
      } else {
        mirna <- annotation2[["name"]]
        mirna_21 <- c(mirna_21, mirna)
      }
    }
  }
  #disconnect DB
  cons <- RMySQL::dbListConnections(RMySQL::MySQL())
  for (con in cons) RMySQL::dbDisconnect(con)
  row.names(data) <- mirna_21
  if (remove_old == TRUE) {
    none_rows <- which(row.names(data) %in% "None")
    if (length(none_rows) != 0) {
      data <- data[-none_rows, ]
    }
  }
  return(data)
}
