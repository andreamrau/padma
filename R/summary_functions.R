#' Summarize results from padma
#'
#' A function to summarize the pathway deviation results from \code{padma}, 
#' using the quantiles of the calculated multi-omic pathway deviation scores.
#'
#' @rdname summary
#' @aliases
#' summary
#' summary-methods
#' summary,padmaResults-method
#'
#' @param object An object of class \code{'padmaResults'}
#' @param ... Additional arguments
#' @author Andrea Rau
#' @seealso \code{\link{padma}}
#' @references
#' Rau, A., Manansala, R., Flister, M. J., Rui, H., Jaffrézic, F., Laloë, D., 
#' and  Auer, P. L. (2019) Individualized multi-omic pathway deviation scores 
#' using multiple factor analysis bioRxiv, https://doi.org/10.1101/827022.
#'
#' @return Summary of the \code{padmaResults} object.
#' @keywords methods
#' @example /inst/examples/padma-package.R
#' @export
setMethod("summary", signature(object = "padmaResults"), 
          function(object, ...) {
    x <- object
    if (!is(x, "padmaResults")) 
        stop(paste0(sQuote("object")), " must be of class ", 
             paste0(dQuote("padmaResults")), 
            sep = "")
    cat("*************************************************\n")
    cat("Pathway multi-omic deviation score quantiles:\n")
    print(quantile(assay(x)$pathway_deviation))
    cat("\nPathway multi-omic deviation scores calculated for:\n", " ", 
        nrow(assay(object)), 
        "individuals\n", " ", nrow(MFA_results(object)$omics), "omics\n", " ", 
        ncol(pathway_gene_deviation(object)), "genes\n")
    cat("*************************************************\n")
})
