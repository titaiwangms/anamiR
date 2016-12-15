#' Find differential expression groups of each genes or miRNA from expression data
#'
#' This function will apply anova ,a statistical methods, for each gene or miRNA (row)
#' to find not only whether expression data of multiple groups differential expressed
#' or not, but also tell specifically two groups from all are differential expression.
#'
#' @return data.frame format with extra columns containing information about
#'    differential expressed groups among all.
#'
#' @seealso \code{\link[stats]{aov}} for fit an analysis of variance model.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'    for input format.
#' @param class string. Choose one features from all rows of phenotype data.
#' @param anova_p_value an numeric value indicating a threshold of p-value from anova
#'    for every genes or miRNAs (rows). Default is 0.05.
#' @param post_hoc post hoc test for anova, including  "scheffe.test", "HSD.test",
#'    "duncan.test".
#' @param post_hoc_p_value an numeric value indicating a threshold of p-value from
#'    post hoc test for every genes or miRNAs (rows). Default is 0.05.
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
#' ## Finding differential miRNA from miRNA expression data with anova
#' aov <- multi_Differ(se = mirna_se, class = "Subtype",
#'    post_hoc = "scheffe.test")
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment assays
#' @import stats
#' @import agricolae
#' @export
multi_Differ <- function(
  se,
  class,
  anova_p_value = 0.05,
  post_hoc = c("scheffe.test",
               "duncan.test",
               "HSD.test"),
  post_hoc_p_value = 0.05
) {

  data <- SummarizedExperiment::assays(se)[[1]]
  pheno_data <- SummarizedExperiment::colData(se)[[class]]

  post_hoc <- match.arg(post_hoc)
  # preprocess
  if (class %in% colnames(pheno_data)) {
    # seperate group
    class.number <- grep(class, colnames(pheno_data))
    group <- as.factor(pheno_data[, class.number])
    group1 <- list()
    # cut data in many pieces
    group_num <- length(levels(group))
    for (i in seq_len(group_num)) {
      group1[[i]] <- data[, which(group == levels(group)[i])]
    }

   # anova
    idx <- c()
    column <- group_num * (group_num - 1) / 2
    count <- matrix(, nrow(data), column)
    for (i in seq_len(nrow(data))) {
      a <- stats::aov(t(data[i, ]) ~ group)
      a_sum <- summary(a)
      if (a_sum[[1]][1, "Pr(>F)"] < anova_p_value) {
        idx <- c(idx, i)
        if (post_hoc %in% "scheffe.test"){
          res <- agricolae::scheffe.test(a,
                                         "group",
                                         group = FALSE,
                                         console = FALSE)
        }
        if (post_hoc %in% "duncan.test"){
          res <- agricolae::duncan.test(a,
                                        "group",
                                        group = FALSE,
                                        console = FALSE)
        }
        if (post_hoc %in% "HSD.test"){
          res <- agricolae::HSD.test(a,
                                     "group",
                                     group = FALSE,
                                     console = FALSE)
        }
        team <- res[["comparison"]]
        for (j in seq_len(nrow(team))) {
          if (team[j, 2] < post_hoc_p_value) {
            count[i, j] <- "TRUE"
          }
        }
      }
    }

    # output
    DE_data <- cbind(data[idx, ], count[idx, ])
    n_samples <- ncol(DE_data)
    same <- row.names(team)
    colnames(DE_data)[(n_samples - (length(same) - 1)):n_samples] <- same
    return(DE_data)
  }
}
