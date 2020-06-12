#' Calculate individualized deviation scores from multi-omic data
#'
#' This is the primary user interface for the \code{padma} package.
#' Generic S4 methods are implemented to calculate individualized pathway
#' deviation scores on the basis of matched, multi-omic data.
#' The supported classes for input are \code{list} and 
#' \code{MultiAssayExperiment}. The output of \code{padma} is an S4 object of 
#' class \code{padmaResults}.
#'
#' @inheritParams padmaRun
#' @param object Matched multi-omic data. May be provided as (1) a 
#' \code{MultiAssayExperiment} or (2) a named \code{list}, with each
#' element corresponding to a \code{data.frame} representing an omic, with 
#' biological entities in rows. Row names should include unique biological 
#' entity IDs (e.g., gene symbols, miRNA names); columns represent 
#' individuals. If more than one biological entity is used, a 
#' \code{gene_map} data.frame providing mappings between IDs and gene
#' names should be provided if the default \code{mirtarbase} is not 
#' sufficient.
#' @param colData (optional) A \code{DataFrame} or \code{data.frame} of 
#' characteristics for all biological units, to be used in creating a 
#' \code{MultiAssayExperiment} from an  \code{object} of class \code{list}
#'
#' @return
#' An S4 object of class \code{padmaResults}, where individualized pathway 
#' deviation scores are stored as the assay data, and the corresponding 
#' {pathway name, full MFA results, number of genes, and names of imputed 
#' or filtered genes} are stored as slots that can be retrieved using
#' the appropriate accessor functions.
#'
#' @aliases
#' padma
#' padma-methods
#' padma,list-method
#' padma,MultiAssayExperiment-method
#'
#' @author Andrea Rau
#' @export
#' @example inst/examples/padma-package.R
#' @keywords methods
#' @rdname padma
#' @docType methods
#' @importFrom methods as is new
setMethod("padma", 
          signature = signature(object = "list"), 
          definition = function(object, colData, gene_map = padma::mirtarbase, 
                                base_ids = NULL, supp_ids = NULL, 
                                pathway_name = "c2_cp_BIOCARTA_D4GDI_PATHWAY", 
                                impute = FALSE, variance_threshold = 1e-04, 
                                full_results = TRUE, verbose = TRUE, ...) {
    omics_data <- object
    arg.user <- list(...)
    
    ## Create MultiAssayExperiment
    if ("character" %in% unique(unlist(
      lapply(object, function(x) apply(x, 2, class))))) 
      stop("One of the data.frames in the provided list contains character data.
  Did you include clinical data in the list? 
  If so, these should instead be omitted or provided via the colData argument.")
    experiments <- object
    omics_data <- suppressMessages(
      MultiAssayExperiment::MultiAssayExperiment(experiments = experiments, 
        colData = colData)
      )
    run <- padmaRun(omics_data = omics_data, gene_map = gene_map, 
                    base_ids = base_ids, 
        supp_ids = supp_ids, pathway_name = pathway_name, impute = impute, 
        variance_threshold = variance_threshold, 
        full_results = full_results, verbose = verbose, ...)
    return(run)
})


############################################################################### 
#' @rdname padma
#' @export
setMethod("padma", 
          signature = signature(object = "MultiAssayExperiment"), 
          definition = function(object, gene_map = padma::mirtarbase, 
                                base_ids = NULL, supp_ids = NULL, 
                                pathway_name = "c2_cp_BIOCARTA_D4GDI_PATHWAY", 
                                impute = FALSE, variance_threshold = 1e-04, 
                                full_results = TRUE, verbose = TRUE, ...) {
    omics_data <- object
    run <- padmaRun(omics_data = omics_data, gene_map = gene_map, 
                    base_ids = base_ids, supp_ids = supp_ids, 
                    pathway_name = pathway_name, impute = impute, 
                    variance_threshold = variance_threshold, 
                    full_results = full_results, verbose = verbose, ...)
    return(run)
})

###############################################################################
#' Accessors for a padmaResults object.
#'
#' @docType methods
#' @rdname padmaHelpers
#' @aliases
#' pathway_name
#' pathway_name,padmaResults-method
#' MFA_results
#' MFA_results,padmaResults-method
#' ngenes
#' ngenes,padmaResults-method
#' imputed_genes
#' imputed_genes,padmaResults-method
#' removed_genes
#' removed_genes,padmaResults-method
#' pathway_gene_deviation
#' pathway_gene_deviation-method
#' show
#' show,padmaResults-method
#' 
#' @param object a \code{padmaResults} object
#' @param ... Additional optional parameters
#' @return Output varies depending on the method. 
#' @author Andrea Rau
#' @export
#' @example inst/examples/padma-package.R
setMethod("pathway_name", signature(object = "padmaResults"), 
          function(object) object@pathway_name)

#' @rdname padmaHelpers
#' @export
setMethod("MFA_results", "padmaResults", function(object) object@MFA_results)

#' @rdname padmaHelpers
#' @export
setMethod("ngenes", "padmaResults", function(object) object@ngenes)

#' @rdname padmaHelpers
#' @export
setMethod("imputed_genes", "padmaResults", 
          function(object) object@imputed_genes)

#' @rdname padmaHelpers
#' @export
setMethod("removed_genes", "padmaResults", 
          function(object) object@removed_genes)

#' @rdname padmaHelpers
#' @export
setMethod("pathway_gene_deviation", "padmaResults", 
          function(object) object@pathway_gene_deviation)

#' @rdname padmaHelpers
#' @export
setMethod("show", 
          signature = "padmaResults", 
          definition = function(object) {
            cat("An object of class ", class(object), "\n", sep = "")
            if (nrow(assay(object))) {
              cat("\nPathway multi-omic deviation scores calculated for:\n", 
                  " ", nrow(assay(object)), "individuals\n", " ", 
                  nrow(MFA_results(object)$omics), "omics\n", " ", 
                  ncol(pathway_gene_deviation(object)), "genes\n\n")
              cat("Use assay(...) to access pathway deviation scores.\n")
    }
    invisible(NULL)
})
