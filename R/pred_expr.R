#' Predicted expression levels for genes in the simulated data.
#'
#' A data set containing predicted expression levels for the 39 genes in the
#' simulated data. We use PTWAS prediction models and a censoring threshold of
#' CPIP > 0.25 to perform expression prediction. This data is included to
#' demonstrate the interface function usage.
#' @format A data frame with 706 rows and 39 variables. Each variable
#' is a gene and each row corresponds to an individual. The row order matches
#' the column order in the included sbams file.
#' @return A data frame with 706 rows and 39 variables.
#' @examples
#' data(pred_expr)
"pred_expr"
