#' Estimate gene priors for causal gene fine-mapping.
#'
#' @param GRCP_vec A vector of colocalization probabilities for the genes of
#' interest
#' @param n_cis A vector containing the number of cis variants for each gene
#' of interest. Order should match that of GRCP_vec.
#' @param p12 A prior enrichment parameter from colocalization analysis.
#' Represents the prior probability that a variant in the region is both a
#' causal GWAS variant and a causal eQTL.
#' @param TWAS_z A vector of TWAS z-scores. Should include all genes
#' genome-wide (likely much larger than GRCP_vec). Order does not matter.
#' @return A vector of gene priors to be used in causal gene fine-mapping.
#' The order of the vector matches the GRCP input.
#' @export
#' @examples
#' data(p1_p12)
#' data(twas_pi1)
#' data(grcp)
#' data(n_cis)
#' isusie_priors(GRCP_vec = grcp, p12 = p1_p12[2],TWAS_pi = twas_pi1,
#' n_cis = n_cis)

isusie_priors <- function(GRCP_vec, p12, TWAS_pi, n_cis){

  p_gene_tmp <- 1-(1-p12)^n_cis

  gene_priors_tmp2 <- 1-(1-GRCP_vec)*(1-p_gene_tmp)

  gene_priors <- gene_priors_tmp2 * TWAS_pi

  return(gene_priors)
}
