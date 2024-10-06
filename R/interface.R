#' Perform probabilistic causal gene fine-mapping with INTERFACE.
#'
#' @param Y An n-vector of GWAS trait measurements for n individuals.
#' @param PE A (n x k) dataframe with predicted expression levels for n
#' individuals and k genes. The order of rows should match Y. Column names
#' should uniquely identify the corresponding gene.
#' @param G A (n x p) dataframe with measured genotypes for n individuals at
#' k variants. The order of rows should match Y. Column names should uniquely
#' identify the corresponding variant.
#' @param p_gene A k-vector of gene priors. We recommend using the
#' interface_priors function to estimate these. The order should match the
#' columns of PE.
#' @param p_snp A scalar variant prior. We recommend using the fastENLOC p1
#' parameter for this quantity.
#' @param null_wt A scalar for the null weight. We recommend estimating this
#' parameter using the interface_null_wt function.
#' @param L An upper bound for the number of causal genes or variants in the
#' region of interest.
#' @return A four-column data frame containing the variable (variant or gene),
#' PIP, posterior effect estimate, and effect estimate posterior standard
#' deviation.
#' @export
#' @importFrom susieR susie
#' @importFrom susieR susie_get_posterior_sd
#' @examples
#' data(sbams)
#' data(pred_expr)
#' data(gene_priors)
#' data(p1_p12)
#' data(null_wt)
#'
#' snp <- sbams$V2[2:nrow(sbams)]
#' G <- t(sbams[2:nrow(sbams),4:ncol(sbams)])
#' colnames(G) <- snp
#' Y = as.numeric(sbams[1,4:709])
#' interface(Y = Y, PE = pred_expr, G = G, p_gene = gene_priors, p_snp = p1_p12[1],null_wt = null_wt)


interface <- function(Y, PE, G, p_gene, p_snp, L = 10, null_wt){

  #Make covariate matrix

  covars <- cbind(PE,G)

  outdf <- data.frame("Variable" = colnames(covars))

  #Perform fine-mapping with SuSiE algorithm

  out <- susie(X = as.matrix(covars),
               y = Y,
               L=10,
               prior_weights = c(p_gene,rep(p_snp,ncol(G)))/(1-null_wt),
               null_weight = null_wt)

  #Extract PIPs, effect estimates, and effect SDs from SuSiE output

  pip <- out$pip

  sds <- as.numeric(susie_get_posterior_sd(out))

  outdf <- cbind(outdf,pip,coef(out)[2:(length(coef(out)) - 1)],
                 sds[1:(length(sds)-1)])

  colnames(outdf) <- c("Variable","PIP","Effect_Estimate","Effect_SD")

  rownames(outdf) <- NULL

  return(outdf)

}

