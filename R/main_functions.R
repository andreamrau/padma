#-----------------------------------------------------------------------
#' Calculate individualized deviation scores from multi-omic data
#'
#' @param omics_data Object of class \code{'MultiAssayExperiment'} containing 
#' omics data from n matched individuals.  
#' @param gene_map (optional) Data frame mapping
#' arbitrary biological entities (e.g. miRNAs) to genes. Contains two columns, 
#' where the first provides the IDs of the entity and
#' the second provides the IDs of the corresponding target gene.
#' By default, the miRNA-gene interactions of type 'Functional MTI' from 
#' miRTarBase are used (see the preloaded \code{'mirtarbase'} data in the 
#' package).
#' @param base_ids (optional) Sample names to be used as reference base data.
#' By default, all samples are used.
#' @param supp_ids (optional) Sample names to be used as supplementary 
#' individuals to be
#' projected onto the analysis based on the individuals identified in 
#' \code{base_ids}. By default, takes the value \code{NULL}, but should not 
#' overlap with \code{base_ids} if provided by the user.
#' @param pathway_name Character of either a KEGG pathway identifier or MSigDB 
#' pathway names (e.g., see the pathway names in the \code{'geneset'} column of 
#' the preloaded \code{msigdb} data in the package), or a vector of gene 
#' symbols.
#' @param impute If \code{TRUE}, impute missing values separately in base and
#' supplementary data using MFA as implemented in the \emph{missMDA} package; 
#' otherwise simple mean imputation is used (default).
#' @param variance_threshold Minimal variance required across samples to retain 
#' a biological entity in the analysis
#' @param full_results If \code{TRUE} (default), include full MFA results in 
#' function output; otherwise, provide concise output to save space.
#' @param verbose If \code{TRUE}, provide verbose output.
#' @param ... Optional additional arguments
#'
#' @return 
#' An S4 object of class \code{padmaResults}, where individualized pathway 
#' deviation scores are stored as the assay data, and the corresponding 
#' {pathway name, full MFA results, number of genes, and names of imputed 
#' or filtered genes} are stored as slots that can be retrieved using
#' the appropriate accessor functions.
#'
#' @export
#' @import MultiAssayExperiment
#' @import SummarizedExperiment
#' @importFrom FactoMineR MFA
#' @importFrom stats var
#' @importFrom utils data
#' @example /inst/examples/padma-package.R
padmaRun <- function(omics_data, gene_map = padma::mirtarbase, base_ids = NULL, 
    supp_ids = NULL, pathway_name = "c2_cp_BIOCARTA_D4GDI_PATHWAY", 
    impute = FALSE, variance_threshold = 1e-04, full_results = TRUE, 
    verbose = TRUE, ...) {
    
    ## Preliminaries --------------------------------------------------------
    if (!is(omics_data, "MultiAssayExperiment")) 
        stop("Expecting a MultiAssayExperiment")
    
    ## Add MSigDB data that are preloaded in the package as msigdb
    MSigDB <- padma::msigdb
    
    ## If base ids are not provided, use all samples (even if user provides
    ## supp_ids)
    if (is.null(base_ids)) {
        base_ids <- rownames(colData(omics_data))
        supp_ids <- NULL
    }
    
    check <- c()
    if (length(supp_ids)) 
        check <- which(base_ids %in% supp_ids)
    if (length(check)) 
        stop("base_ids and supp_ids should be mutually exclusive.")
    
    ## Identify pathway genes
    ## -----------------------------------------------------------------
    if (length(pathway_name) == 1) {
        ## KEGG IDs
        if (length(grep("hsa", pathway_name))) {
            if (!requireNamespace("KEGGREST")) 
stop("Please load KEGGREST to automatically search for KEGG pathway gene IDs.")
            kg <- KEGGREST::keggGet(pathway_name)[[1]]
            Pathway_Name <- paste0(pathway_name, ": ", kg$NAME)
            pathway_genes <- data.frame(gene = kg$GENE)
            if (!nrow(pathway_genes)) {
                RESULTS <- NULL
                return(RESULTS)
            } else {
                pathway_genes <- unique(unlist(lapply(
                    strsplit(pathway_genes[grep(";", 
                  pathway_genes[, 1]), 1], split = "; "), function(x) x[1])))
            }
            ## C2 built-in MSigDB pathways
        } else if ((length(grep("c2_cp_", pathway_name)))) {
            Pathway_Name <- pathway_name
            pathway_genes <- unlist(strsplit(
                MSigDB[which(MSigDB$geneset == pathway_name), 
                2], split = ", "))
        } else {
            stop("No other built-in pathways currently supported.")
        }
        ## Custom pathways
    } else {
        ## pathway_name consists of a vector of gene names
        Pathway_Name <- "Custom"
        pathway_genes <- pathway_name
    }
    
    ## Subset data -------------------------------------------------------
    expanded_pathway_genes <- pathway_genes
    if (!is.null(gene_map)) {
        gene_map_subset <- gene_map[which(gene_map[, 2] %in% pathway_genes), , 
            drop = FALSE]
        gene_map_subset[, 1] <- tolower(gene_map_subset[, 1])
        entity_subset <- gene_map_subset[, 1]
        expanded_pathway_genes <- c(pathway_genes, entity_subset)
    }
    base_data <- suppressMessages(omics_data[expanded_pathway_genes, base_ids, 
        , drop = FALSE])
    supp_data <- suppressMessages(omics_data[expanded_pathway_genes, supp_ids, 
        , drop = FALSE])
    
    ## Add gene names to other entities (and duplicate) as needed
    if (!is.null(gene_map)) {
        for (j in names(base_data)) {
            ids <- rownames(experiments(base_data)[[j]])
            if (sum(ids %in% gene_map_subset[, 1])) {
                ids_tmp <- unique(
                    gene_map_subset[which(gene_map_subset[, 1] %in% 
                  ids), , drop = FALSE])
                exp_tmp <- matrix(
                  experiments(base_data)[[j]][match(ids_tmp[,1], 
                       rownames(experiments(base_data)[[j]])), ],
                  ncol = ncol(experiments(base_data)[[j]]))
                colnames(exp_tmp) <- colnames(experiments(base_data)[[j]])
                experiments(base_data)[[j]] <- exp_tmp
                rownames(experiments(base_data)[[j]]) <- paste0(ids_tmp[, 1], 
                  "_", ids_tmp[, 2])
            }
        }
    }
    
    ## Only keep samples that are represented by all assays
    base_data <- suppressMessages(intersectColumns(base_data))
    ## OCt 3, 2022 fix: intersectColumns does not work for empty mae
    if(max(unlist(lapply(assays(supp_data), ncol))) > 0) {
      supp_data <- suppressMessages(intersectColumns(supp_data))
    }

    ## Return empty results if no genes
    if(max(unlist(lapply(experiments(base_data), nrow))) == 0) {
      message("No matching genes found.")
      RESULTS <- NULL
      return(RESULTS)
    }
    
    ## Filter / impute data for entities as needed ----------------------
    nc_base <- ncol(base_data[[1]])
    nc_supp <- ncol(supp_data[[1]])
    removed_genes <- imputed_genes <- vector("list", length(names(base_data)))
    names(removed_genes) <- names(imputed_genes) <- names(base_data)
    ## If there are only base individuals, remove genes that have 0 var or
    ## all NAs
    if (!nrow(sampleMap(supp_data))) {
        for (j in names(base_data)) {
            remove_index <- unique(c(which(apply(base_data[[j]], 1, var) < 
                                               variance_threshold), 
                which(rowSums(is.na(base_data[[j]])) == nc_base)))
            if (length(remove_index)) {
                removed_genes[[j]] <- unlist(rownames(
                    base_data[[j]])[remove_index])
                experiments(base_data)[[j]] <- 
                    experiments(base_data)[[j]][-remove_index, 
                  ]
                if (verbose) {
                  cat("The following are removed from", j, 
                      "in base due to small variance or all NAs:\n")
                  print(removed_genes[[j]])
                }
            }
        }
    } else {
        ## If base + supp individuals, remove genes with small variance / NA
        ## in both. For genes with small var in base alone, impute genes with 0
        ## variability by reshuffling the minimally variant (> 10e-5) value
        for (j in names(base_data)) {
            ## Remove elements with very small variance or all NA's in both
            remove_base <- unique(c(which(apply(base_data[[j]], 1, var, 
                                                na.rm = TRUE) < 
                variance_threshold), which(rowSums(is.na(base_data[[j]])) == 
                                               nc_base)))
            remove_supp <- unique(c(which(apply(supp_data[[j]], 1, var, 
                                                na.rm = TRUE) < 
                variance_threshold), which(rowSums(is.na(supp_data[[j]])) == 
                                               nc_supp)))
            remove_index <- intersect(remove_base, remove_supp)
            if (length(remove_index)) {
                removed_genes[[j]] <- unlist(rownames(
                    base_data[[j]])[remove_index])
                experiments(base_data)[[j]] <- 
                    experiments(base_data)[[j]][-remove_index, 
                  ]
                experiments(supp_data)[[j]] <- 
                    experiments(supp_data)[[j]][-remove_index, 
                  ]
                if (verbose) {
                  cat("The following were removed from", j, 
                      "in base and supplementary due to small variance:\n")
                  print(removed_genes[[j]])
                }
            }
            
            ## Now identify elements with very small variance or all NA's 
            ## in base alone to impute
            impute_index <- unique(c(which(apply(base_data[[j]], 1, var) < 
                                               variance_threshold), 
                which(rowSums(is.na(base_data[[j]])) == nc_base)))
            if (length(impute_index)) {
                imputed_genes[[j]] <- unlist(rownames(
                    base_data[[j]])[impute_index])
                
                if (verbose) {
                  cat("The following were imputed for", j, 
                      "in base due to small variance or all NAs:\n")
                  print(imputed_genes[[j]])
                }
                wmv <- suppressWarnings(whichminpositivevar(base_data[[j]]))
                if (!is.na(wmv)) {
                  choose_impute <- unlist(base_data[[j]][wmv, ])
                } else {
                  choose_impute <- c(1, rep(0, nc_base - 1))
                }
                for (jj in seq_len(length(impute_index))) {
                  base_data[[j]][impute_index[jj], ] <- sample(choose_impute)
                }
            }
            ## Fill in remaining all NA's in supp with 0's
            supp0_index <- which(rowSums(is.na(supp_data[[j]])) == nc_supp)
            if (length(supp0_index) & nc_supp > 0) {
                supp_data[[j]][supp0_index, ] <- 0
            }
        }
    }
    
    ## Formatting data: use wideFormat ----------------------------------------
    gene_tables_base <- wideFormat(base_data, check.names = FALSE)
    colnames(gene_tables_base) <- switch_names(colnames(gene_tables_base))
    rownames(gene_tables_base) <- gene_tables_base[, 1]
    gene_tables_base <- gene_tables_base[, -1]
    ## Reorder to group by genes
    o <- order(colnames(gene_tables_base))
    gene_tables_base <- gene_tables_base[, o]
    if (nrow(colData(supp_data))) {
        gene_tables_supp <- wideFormat(supp_data, check.names = FALSE)
        colnames(gene_tables_supp) <- switch_names(colnames(gene_tables_supp))
        rownames(gene_tables_supp) <- gene_tables_supp[, 1]
        gene_tables_supp <- gene_tables_supp[, -1]
        ## Reorder to group by genes
        o <- order(colnames(gene_tables_supp))
        gene_tables_supp <- gene_tables_supp[, o]
    } else {
        gene_tables_supp <- gene_tables_base[-seq_len(nrow(gene_tables_base)), 
            ]
    }
    
    ## Remove elements with more missing elements than half the number of base
    ## individuals in either group
    remove_index <- unique(c(which(colSums(is.na(gene_tables_base)) > 0.5 * 
                                       nrow(base_data$rnaseq))))
    if (length(remove_index)) {
        gene_tables_supp <- gene_tables_supp[, -remove_index]
        gene_tables_base <- gene_tables_base[, -remove_index]
    }
    
    ## Data imputation ------------------------------------------------------
    ## Only use missMDA for data imputation if specified by the user !!!  
    ## Impute missing data for base and supp independently with missMDA 
    ## (use type = 's' to scale data)
    lgr_tab <- table(unlist(lapply(strsplit(colnames(gene_tables_base), 
                                            split = "_", 
        fixed = TRUE), function(x) x[1])))
    orig <- unique(unlist(lapply(strsplit(colnames(gene_tables_base), 
                                          split = "_", 
        fixed = TRUE), function(x) x[1])))
    lgr_tab <- lgr_tab[match(orig, names(lgr_tab))]
    lgr <- as.numeric(lgr_tab)
    names(lgr) <- names(lgr_tab)
    if (sum(is.na(gene_tables_base)) & impute) {
        if (!requireNamespace("missMDA")) 
      stop("Please load missMDA to impute missing data using using inputeMFA.")
        gene_tables_base <- missMDA::imputeMFA(gene_tables_base, group = lgr, 
            ncp = 2, type = rep("s", 
                                length(names(gene_tables_base))))$completeObs
    }
    if (sum(is.na(gene_tables_supp)) & impute) {
        if (!requireNamespace("missMDA")) 
     stop("Please load missMDA to impute missing data using using inputeMFA.")
        gene_tables_supp <- missMDA::imputeMFA(gene_tables_supp, group = lgr, 
            ncp = 2, type = rep("s", 
                                length(names(gene_tables_supp))))$completeObs
    }
    
    gene_tables_base_s <- scale(gene_tables_base, center = TRUE, scale = TRUE)
    if (nrow(gene_tables_supp)) {
        gene_tables_supp <- as(gene_tables_supp, "data.frame")
        gene_tables_supp_s <- 
            t((t(gene_tables_supp) - 
                   attributes(gene_tables_base_s)$`scaled:center`) / 
                  attributes(gene_tables_base_s)$`scaled:scale`)
    } else {
        gene_tables_supp_s <- gene_tables_supp
    }
    
    ## Combine all observations
    if (nrow(gene_tables_supp_s)) {
        gene_tables <- rbind(gene_tables_base_s, gene_tables_supp_s)
    } else {
        gene_tables <- gene_tables_base_s
    }
    
    ## Run MFA ----------------------------------------------------------------
    ## c = pre-scaled variables
    group <- c(rep("Base", nrow(gene_tables_base_s)), 
               rep("Supp", nrow(gene_tables_supp_s)))
    if (!nrow(gene_tables_supp_s)) {
        ind.sup <- NULL
    } else {
        ind.sup <- (nrow(gene_tables_base) + 1):(nrow(gene_tables_base) + 
                                                     nrow(gene_tables_supp))
    }
    total_MFA <- MFA(gene_tables, group = as.vector(lgr), 
                     type = as.vector(rep("c", 
        length(lgr))), ind.sup = ind.sup, ncp = ncol(gene_tables), 
        graph = FALSE, 
        name.group = names(lgr))
    
    ## Calculate downstream scores---------------------------------------------
    
    ## Calculate pathway deregulation score
    ps <- pathway_scores(total_MFA, ngenes = length(lgr))
    ps_df <- data.frame(primary = names(ps), group = group, 
                        pathway_deviation = ps, 
        row.names = NULL, stringsAsFactors = FALSE)
    
    ## Calculate the gene-specific pathway deregulation scores
    f <- total_MFA$ind$coord
    fk <- total_MFA$ind$coord.partiel
    gene_names <- unique(unlist(lapply(strsplit(colnames(gene_tables), "_", 
                                                fixed = TRUE), 
        function(x) x[1])))
    fk_list <- vector("list", length(gene_names))
    names(fk_list) <- gene_names
    for (g in gene_names) {
        fk_list[[g]] <- fk[grep(paste0(".", g, "$"), rownames(fk)), ]
    }
    df_final <- matrix(NA, nrow = nrow(f), ncol = length(gene_names))
    rownames(df_final) <- rownames(f)
    colnames(df_final) <- gene_names
    for (g in gene_names) {
        ## Renormalized by distance scores
        df_final[, g] <- rowSums((f * (fk_list[[g]] - f)))/sqrt(rowSums(f^2))
    }
    gs_df <- df_final
    
    ## Omics: % contribution to each axis of the MFA
    omics_contrib_MFA <- NULL
    omics_groups <- unlist(lapply(strsplit(rownames(
        total_MFA$quanti.var$contrib), 
        split = "_", fixed = TRUE), function(x) paste0(x[-1], collapse = "_")))
    
    if (!is.na(omics_groups[1])) {
        if (length(grep("_", omics_groups))) {
            omics_groups[grep("_", omics_groups)] <- 
                unlist(lapply(strsplit(omics_groups[grep("_", 
                omics_groups)], split = "_", fixed = TRUE), function(x) x[2]))
        }
        omics_contrib_MFA <- round(rowsum(total_MFA$quanti.var$contrib, 
                                          group = omics_groups), 2)
    }
    
    if (!full_results) {
        if (!is.na(omics_groups[1])) {
            omics_contrib_MFA_summary <- data.frame(
                PC_10 = round(rowSums(omics_contrib_MFA[, 
                  seq_len(min(10, ncol(omics_contrib_MFA)))]) / 
                      min(10, ncol(omics_contrib_MFA)), 2), 
                PC_all = round(rowSums(
                    omics_contrib_MFA[, seq_len(ncol(omics_contrib_MFA))]) / 
                        ncol(omics_contrib_MFA), 2))
        } else {
            omics_contrib_MFA_summary <- data.frame(PC_10 = 1, PC_all = 1)
            row.names(omics_contrib_MFA_summary) <- "single-omic"
        }
    }
    
    ## Format results ---------------------------------------------------------
    pathway_deviation <- SummarizedExperiment(ps_df)
    if (full_results) {
        ## Base individuals: perc contributions to each PC Genes: perc 
        ## contributions to each axis of the MFA Genes: Lg coefficients to show 
        ## pairwise correlations among genes Omics: perc contribution to each 
        ## axis of the MFA Full MFA results and gene tables used
        
        MFA_results <- list(eig = total_MFA$eig, 
                            ind_contrib_MFA = round(total_MFA$ind$contrib, 
            2), gene_contrib_MFA = round(total_MFA$group$contrib, 2), 
            gene_Lg_MFA = round(total_MFA$group$Lg, 
            4), omics_contrib_MFA = omics_contrib_MFA, total_MFA = total_MFA, 
            gene_tables = gene_tables)
    } else {
        ctrb <- total_MFA$ind$contrib
        ind <- min(10, ncol(ctrb))
        grp <- total_MFA$group$contrib
        ind2 <- min(10, ncol(grp))
        ## Base individuals: 
        ##   perc contributions to each first 10 PCs and all PCs Genes:
        ## perc contributions to each first 10 PCs and all PCs Omics: 
        ##   perc contribution to each first 10 PCs and all PCs
        MFA_results <- 
            list(eig = total_MFA$eig, 
                 ind_contrib_MFA_summary = data.frame(
                     PC_10 = round(rowSums(ctrb[, seq_len(ind)])/ind, 2), 
                     PC_all = round(rowSums(
                         ctrb[, seq_len(ncol(ctrb))])/ncol(ctrb), 2)), 
                 gene_contrib_MFA_summary = data.frame(
                     PC_10 = round(rowSums(grp[,seq_len(ind2)])/ind2, 2), 
                     PC_all = round(rowSums(grp[, seq_len(ncol(grp))]) / 
                                        ncol(grp), 2)), 
                 omics_contrib_MFA_summary = omics_contrib_MFA_summary)
    }
    
    RESULTS <- padmaResults(as(pathway_deviation, 
                               "RangedSummarizedExperiment"), 
        pathway_name = Pathway_Name, 
        pathway_gene_deviation = DataFrame(gs_df), 
        ngenes = length(lgr), imputed_genes = imputed_genes, 
        removed_genes = removed_genes, 
        MFA_results = MFA_results)
    return(RESULTS)
}

#-----------------------------------------------------------------------

## DO NOT EXPORT
pathway_scores <- function(MFA_results, ngenes) {
    return(deviation = sqrt(rowSums(
        rbind(MFA_results$ind$coord, MFA_results$ind.sup$coord)^2)))/ngenes
}

## DO NOT EXPORT
whichminpositivevar <- function(xx, threshold = 1e-04) {
    xxx <- apply(xx, 1, var, na.rm = TRUE)
    return(which(xxx == min(xxx[xxx > threshold], na.rm = TRUE))[1])
}

## DO NOT EXPORT
switch_names <- function(name_vector) {
    tmp <- unlist(lapply(lapply(strsplit(name_vector, split = "_", 
                                         fixed = TRUE), 
        rev), function(x) paste0(x, collapse = "_")))
    return(tmp)
}
