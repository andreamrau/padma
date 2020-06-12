#' @rdname padma
#' @export
setGeneric("padma", function(object, ...) standardGeneric("padma"))

#' @rdname padmaHelpers
#' @export
setGeneric("pathway_name", 
           function(object, ...) standardGeneric("pathway_name"))

#' @rdname padmaHelpers
#' @export
setGeneric("pathway_gene_deviation", 
           function(object, ...) standardGeneric("pathway_gene_deviation"))

#' @rdname padmaHelpers
#' @export
setGeneric("MFA_results", 
           function(object, ...) standardGeneric("MFA_results"))

#' @rdname padmaHelpers
#' @export
setGeneric("ngenes", function(object, ...) standardGeneric("ngenes"))

#' @rdname padmaHelpers
#' @export
setGeneric("imputed_genes", 
           function(object, ...) standardGeneric("imputed_genes"))

#' @rdname padmaHelpers
#' @export
setGeneric("removed_genes", 
           function(object, ...) standardGeneric("removed_genes"))
