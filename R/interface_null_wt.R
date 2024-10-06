#' Estimate null weight for causal gene fine-mapping
#'
#' @param p1 A prior enrichment parameter from colocalization analysis.
#' Represents the prior probability that a variant in the region is a
#' causal GWAS variant but not a causal eQTL.
#' @param p12 A prior enrichment parameter from colocalization analysis.
#' Represents the prior probability that a variant in the region is both a
#' causal GWAS variant and a causal eQTL.
#' @param n_region_snps The number of variants in the region of interest for
#' fine-mapping.
#' @return The null weight (scalar).
#' @export
#' @examples
#' data(p1_p12)
#' data(sbams)
#'
#' interface_null_wt(p1 = p1_p12[1], p12 = p1_p12[2],n_region_snps=(nrow(sbams)-1))

interface_null_wt <- function(p1, p12, n_region_snps){

  p_gwas = p1 + p12

  null_wt = (1-p_gwas)^n_region_snps

  return(null_wt)
}
