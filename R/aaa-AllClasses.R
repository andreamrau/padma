#' @rdname padmaResults
#' @export
#' @import SummarizedExperiment
setClass("padmaResults", contains = "RangedSummarizedExperiment", 
         representation = 
             representation(pathway_name = "character", 
                            pathway_gene_deviation = "DataFrame", 
                            MFA_results = "list", 
                            ngenes = "integer", 
                            imputed_genes = "list", 
                            removed_genes = "list"))


#' padmaResults object and constructor
#'
#' \code{padmaResults} is a subclass of \code{RangedSummarizedExperiment},
#' used to store the individualized pathway deviation scores as well as some 
#' additional information useful about the pathway name (\code{pathway_name}), 
#' the gene-level contributions to each deviation score 
#' (\code{pathway_gene_deviation}), a full set of 
#' outputs related to the MFA (\code{MFA_results}, and the number
#' of genes used in the analysis as well as the names of those for which data 
#' imputation or filtering was required (\code{ngenes}, \code{imputed_genes}, 
#' and \code{removed_genes}, respectively).
#'
#' This constructor function would not typically be used by 'end users'.
#' This simple class extends the \code{RangedSummarizedExperiment} class of the
#' SummarizedExperiment package
#' to allow other packages to write methods for results
#' objects from the padma package. It is used by \code{\link{padmaRun}}
#' to wrap up the results table.
#'
#' @param SummarizedExperiment a \code{RangedSummarizedExperiment} of 
#' \code{padma} results
#' @param pathway_name The name of the pathway, if applicable
#' @param pathway_gene_deviation Per-gene contributions to each 
#' individualized pathway deviation score 
#' @param MFA_results List of all detailed results from the MFA
#' @param ngenes Number of genes used in the pathway deviation score calculation
#' @param imputed_genes Names of genes, per omic, for which data imputation was 
#' used to replace missing values
#' @param removed_genes Names of genes, per omic, which were filtered from the 
#' analysis due to low variation
#' @return a padmaResults object
#' @docType class
#' @rdname padmaResults
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom S4Vectors DataFrame
padmaResults <- function(SummarizedExperiment, 
                         pathway_name = NULL, 
                         pathway_gene_deviation = NULL, 
                         MFA_results = NULL, 
                         ngenes = NULL, 
                         imputed_genes = NULL, 
                         removed_genes = NULL) {
    se <- SummarizedExperiment
    if (!is(se, "RangedSummarizedExperiment")) {
     stop("'SummarizedExperiment' must be a RangedSummarizedExperiment object")
    }
    if (is.null(pathway_name)) 
        pathway_name <- ""
    if (is.null(pathway_gene_deviation)) 
        pathway_gene_deviation <- DataFrame(matrix(0, nrow = 0, ncol = 0))
    if (is.null(MFA_results)) 
        MFA_results <- list()
    if (is.null(ngenes)) 
        ngenes <- integer()
    if (is.null(imputed_genes)) 
        imputed_genes <- list()
    if (is.null(removed_genes)) 
        removed_genes <- list()
    
    object <- new("padmaResults", se, pathway_name = pathway_name, 
                  pathway_gene_deviation = pathway_gene_deviation, 
                  MFA_results = MFA_results, ngenes = ngenes, 
                  imputed_genes = imputed_genes, removed_genes = removed_genes)
    return(object)
}
