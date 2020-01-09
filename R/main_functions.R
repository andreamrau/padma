
## TODO:
## - Generalize to arbitrary names of list
## - Generalize to non-TCGA sample IDs
## - Fix "no visible binding for global variable" errors due to tidyverse

#-----------------------------------------------------------------------
#' Calculate individualized deviation scores from multi-omic data
#'
#' @param omics_data Named list (including "clinical" and "rnaseq", and optionally
#' additionaly including one or more of "cna", "methyl", and/or "mirna"),
#' where each element of the list correspnds to a data.frame with
#' omics data from n matched individuals. The data.frame for clinical
#' data ("clinical") has dimension n x (q + 1) for q clinical variables, where an additional
#' column named "bcr_patient_barcode" contains unique TCGA sample identifiers.
#' The data.frames for
#' RNA-seq ("rnaseq"), copy number alterations ("cna"), and methylation ("methyl") have
#' dimension n x (p + 1) for p genes, where the first column contains gene names. The
#' data.frame for microRNA-seq ("mirna") has dimension n x (p + 2) for p genes, where the
#' first column contains the miRNA ID and the second contains the gene name of predicted
#' targets. Row names of all data.frames should contain the same unique sample identifiers
#' as those included in the "bcr_patient_barcode" variable of the "clinical" data.frame.
#' Sample orders (rows of "clinical" and columns of the others) should match across all
#' data.frames
#' @param base_ids (optional) Sample indices to be used as reference base data.
#' By default, all samples are used.
#' @param supp_ids (optional) Sample indices to be used as supplementary individuals to be
#' projected onto the analysis based on the individuals identified in \code{base_ids}. By
#' default, takes the value \code{NULL}, but should not overlap with \code{base_ids} if
#' provided by the user.
#' @param apply_log (optional) Data types (corresponding to the names of the \code{omics_data}
#' list) to which a log + 1 transformation should be applied.
#' @param pathway_name Either a character of a KEGG pathway identifier or MSigDB pathway names
#' (e.g., see the pathway names in the \code{"geneset"} column of the preloaded
#' \code{msigdb} data in the package), or a vector of gene symbols.
#' @param impute_MFA_missMDA If \code{TRUE}, impute missing values separately in base and
#' supplementary data using MFA as implemented in the \emph{missMDA} package; otherwise mean
#' imputation is used (default).
#' @param mirna_targets (optional) Data.frame with 2 columns: column \code{"miRNA"}
#' should provide a miRNA name (e.g., \code{"hsa-miR-20a-5p"}) and column \code{"Target Gene"}
#' should provide the symbol of a predicted gene target for the corresponding miRNA.
#' Each predicted miRNA-target pair should have its own row in the
#' data.frame. By default, the miRNA-gene interactions of type "Functional MTI" from miRTarBase
#' are used (see the preloaded \code{"mirtarbase"} data in the package).
#' @param full_results If \code{TRUE} (default), include full MFA results in function output;
#' otherwise, provide concise output to save space.
#'
#' @return Output of class \code{"padma"} containing a list of several elements.
#'
#' @export
#' @import tidyr dplyr
#' @importFrom KEGGREST keggGet
#' @importFrom missMDA imputeMFA
#' @importFrom FactoMineR MFA
#' @importFrom stats quantile var weighted.mean
#' @importFrom utils data tail
#' @example /inst/examples/padma-package.R
padma <- function(omics_data,
                  base_ids = NULL,
                  supp_ids = NULL,
                  apply_log = NULL,
                  pathway_name = "c2_cp_BIOCARTA_D4GDI_PATHWAY",
                  impute_MFA_missMDA = FALSE,
                  mirna_targets = NULL,
                  full_results = TRUE) {

  ## Add MSigDB data that are preloaded in the package as misgdb
  MSigDB <- padma::msigdb
  ## Add miRTarBase mirna_targets preloaded as mirtarbase if mirna_targets not provided
  miRTarBase <- padma::mirtarbase
  if(!is.null(mirna_targets)) {
    miRTarBase <- mirna_targets
  }
  ## If base ids are not provided, use all samples (even if user provides supp_ids)
  if(is.null(base_ids)) {
    base_ids <- seq_len(nrow(omics_data$clinical))
    supp_ids <- NULL
  }

  check <- c()
  if(length(supp_ids)) check <- which(base_ids %in% supp_ids)
  if(length(check)) stop("base_ids and supp_ids should be mutually exclusive.")
  if(is.null(omics_data$clinical$bcr_patient_barcode)) stop("clinical data should include bcr_patient_barcode variable.")
  if(!"clinical" %in% names(omics_data)) stop("clinical data must be included.")
  if(!"rnaseq" %in% names(omics_data)) stop("rnaseq data must be included.")

  ## Select pathway genes ---------------------------------------------------------------------------
  if(length(pathway_name) == 1) {
    if(length(grep("hsa", pathway_name))) {        ## KEGG IDs
      kg <- keggGet(pathway_name)[[1]]
      Pathway_Name <- paste0(pathway_name, ": ", kg$NAME)
      pathway_genes <- data.frame(gene = kg$GENE)
      if(nrow(pathway_genes) == 0) {
        return(list(Pathway_Name=Pathway_Name,
                    pathway_deregulation = NULL,
                    gene_deregulation = NULL,
                    omics_deregulation = NULL,
                    eig = NULL,
                    coordinates_MFA = NULL))
      } else {
        pathway_genes <- pathway_genes %>%
          slice(grep(";", gene)) %>%
          separate(gene, into = c("gene","details"), sep="; ") %>%
          dplyr::select(gene) %>%
          unique() %>% unlist()
      }
    } else if((length(grep("c2_", pathway_name)))) { ## C2 MSigDB pathways
      Pathway_Name <- pathway_name
      pathway_genes <- MSigDB %>%
        dplyr::filter(geneset == pathway_name) %>%
        dplyr::select(symbol) %>% unlist() %>% strsplit(., split = ", ") %>% unlist()
    } else {
      stop("No other pathways currently supported.")
    }
  } else {
    ## pathway_name consists of a vector of gene names
    Pathway_Name <- "Custom"
    pathway_genes <- pathway_name
  }

  ## Subset data ------------------------------------
  ps <- pathway_subset(x = omics_subset(omics_data, index = base_ids, apply_log = apply_log),
                       y = omics_subset(omics_data, index = supp_ids, apply_log = apply_log),
                       subset = pathway_genes, miRTarBase = miRTarBase)
  base_data <- ps$x
  supp_data <- ps$y

  ## Stop if there are fewer than 3 observed genes in the pathway
  if(nrow(base_data$rnaseq) < 3) {
    return(list(Pathway_Name=Pathway_Name,
                pathway_deregulation = NULL,
                gene_deregulation = NULL,
                omics_deregulation = NULL,
                eig = NULL,
                coordinates_MFA = NULL))
  }

  removed_genes <- imputed_genes <- vector("list", 4)
  names(removed_genes) <- names(imputed_genes) <- c("rnaseq", "methyl", "mirna", "cna")
  ## If there are only base individuals, remove genes that have 0 variability or all NAs
  if(ncol(supp_data$rnaseq) == 1) {  ## First column of RNA-seq is the gene name
    for(j in names(base_data)) {
      if(!is.null(base_data[[j]]) & j != "clinical") {
        if(j == "mirna") {
          rem_col <- c(1,2)
        } else {
          rem_col <- 1
        }
        remove_index <- unique(c(which(apply(base_data[[j]][,-rem_col], 1, var) < 10e-5),
                                 which(rowSums(is.na(base_data[[j]][,-rem_col])) == ncol(base_data[[j]][,-rem_col]))))
        if(length(remove_index) > 0) {
          cat("The following are removed from", j, "in base due to small variance or all NAs:\n")
          print(base_data[[j]][remove_index,1])
          unlist(removed_genes[[j]] <- base_data[[j]][remove_index,1])
          base_data[[j]] <- base_data[[j]][-remove_index,]
        }
      }
    }
  }
  ## If there are base + supp individuals, remove genes with small variance in both base and supp.
  ## For genes with small variance in base alone, "impute" genes with 0 variability by randomly
  ##   reshuffling the minimally variant (> 10e-5) value
  if(ncol(supp_data$rnaseq) > 1) {
    for(j in names(base_data)) {
      if(!is.null(base_data[[j]]) & j != "clinical") {
        if(j == "mirna") {
          rem_col <- c(1,2)
        } else {
          rem_col <- 1
        }
        ## Remove elements with very small variance or all NA's in both
        remove_base <- unique(c(which(apply(base_data[[j]][,-rem_col], 1, var, na.rm=TRUE) < 10e-5),
                                which(rowSums(is.na(base_data[[j]][,-rem_col])) == ncol(base_data[[j]][,-rem_col]))))
        remove_supp <- unique(c(which(apply(supp_data[[j]][,-rem_col], 1, var, na.rm=TRUE) < 10e-5),
                                which(rowSums(is.na(supp_data[[j]][,-rem_col])) == ncol(supp_data[[j]][,-rem_col]))))
        remove_index <- intersect(remove_base, remove_supp)
        if(length(remove_index)) {
          cat("The following were removed from", j, "in base and supplementary due to small variance:\n")
          print(unlist(base_data[[j]][remove_index,1]))
          removed_genes[[j]] <- unlist(base_data[[j]][remove_index,1])
          base_data[[j]] <- base_data[[j]][-remove_index,]
          supp_data[[j]] <- supp_data[[j]][-remove_index,]
        }
        ## Now identify elements with very small variance or all NA's in base alone to impute
        impute_index <- unique(c(which(apply(base_data[[j]][,-rem_col], 1, var) < 10e-5),
                                 which(rowSums(is.na(base_data[[j]][,-rem_col])) == ncol(base_data[[j]][,-rem_col]))))
        if(length(impute_index)) {
          cat("The following were imputed for", j, "in base due to small variance or all NAs:\n")
          print(unlist(base_data[[j]][impute_index,1]))
          imputed_genes[[j]] <- unlist(base_data[[j]][impute_index,1])

          wmv <- suppressWarnings(whichminpositivevar(base_data[[j]][,-rem_col]))
          if(!is.na(wmv)) {
            choose_impute <- unlist(base_data[[j]][wmv,-rem_col])
          } else {
            choose_impute <- c(1,rep(0, ncol(base_data[[j]][,-rem_col])-1))
          }
          for(jj in 1:length(impute_index)) {
            if(j %in% c("rnaseq", "mirna")) {
              base_data[[j]][impute_index[jj],-rem_col] <- sample(choose_impute)
            } else {
              base_data[[j]][impute_index[jj],-rem_col] <- scale(sample(choose_impute),
                                                                 center=TRUE, scale=FALSE)
            }
          }
        }
        ## Fill in remaining all NA's in supp with 0's
        supp0_index <- which(rowSums(is.na(supp_data[[j]][,-rem_col])) == ncol(supp_data[[j]][,-rem_col]))
        if(length(supp0_index)) {
          supp_data[[j]][supp0_index,-rem_col] <- 0
        }
      }
    }
  }

  ## Formatting data
  gene_tables_supp <- gene_tables_base <- vector("list", nrow(base_data$rnaseq))
  names(gene_tables_supp) <- names(gene_tables_base) <- base_data$rnaseq$gene
  for(i in names(gene_tables_base)) {
    gene_tables_base[[i]] <- data.frame(rnaseq = dplyr::filter(base_data$rnaseq, gene == i) %>%
                                          dplyr::select(-gene) %>% unlist(),
                                        check.names=FALSE, stringsAsFactors = FALSE)
    gene_tables_supp[[i]] <- data.frame(rnaseq = dplyr::filter(supp_data$rnaseq, gene == i) %>%
                                          dplyr::select(-gene) %>% unlist(),
                                        check.names=FALSE, stringsAsFactors = FALSE)
    ## Add in methylation data if available
    if(!is.null(base_data$methyl)) {
      if(nrow(dplyr::filter(base_data$methyl, gene == i)) > 0) { ## AR: check if methyl is available for that gene

        gene_tables_base[[i]] <- data.frame(gene_tables_base[[i]], methyl=dplyr::filter(base_data$methyl, gene == i) %>%
                                              dplyr::select(-gene) %>% unlist(),
                                            check.names=FALSE, stringsAsFactors = FALSE)
        gene_tables_supp[[i]] <- data.frame(gene_tables_supp[[i]], methyl=dplyr::filter(supp_data$methyl, gene == i) %>%
                                              dplyr::select(-gene) %>% unlist(),
                                            check.names=FALSE, stringsAsFactors = FALSE)
      }
    }
    ## Add in miRNA data if available
    if(!is.null(base_data$mirna)) {
      if(nrow(dplyr::filter(base_data$mirna, `Target Gene` == i)) > 0) { ## AR: check if miR available for that gene
        tmp <- dplyr::filter(base_data$mirna, `Target Gene` == i) %>%
          dplyr::select(-`Target Gene`)
        mirs <- as.character(tmp$miRNA_lc)
        tmp <- tmp %>% dplyr::select(-miRNA_lc) %>% t()
        colnames(tmp) <- mirs
        gene_tables_base[[i]] <- data.frame(gene_tables_base[[i]], tmp,
                                            check.names=FALSE, stringsAsFactors = FALSE)
        tmp <- dplyr::filter(supp_data$mirna, `Target Gene` == i) %>%
          dplyr::select(-`Target Gene`)
        mirs <- as.character(tmp$miRNA_lc)
        tmp <- tmp %>% dplyr::select(-miRNA_lc) %>% t()
        colnames(tmp) <- mirs
        gene_tables_supp[[i]] <- data.frame(gene_tables_supp[[i]], tmp,
                                            check.names=FALSE, stringsAsFactors = FALSE)
      }
    }
    ## Add in CNA data if available
    if(!is.null(base_data$cna)) {
      if(nrow(dplyr::filter(base_data$cna, gene == i)) > 0) { ## AR: check if CNA is available for that gene
        gene_tables_base[[i]] <- data.frame(gene_tables_base[[i]], cna=dplyr::filter(base_data$cna, gene == i) %>%
                                              dplyr::select(-gene) %>% unlist(),
                                            check.names=FALSE, stringsAsFactors = FALSE)
        gene_tables_supp[[i]] <- data.frame(gene_tables_supp[[i]], cna=dplyr::filter(supp_data$cna, gene == i) %>%
                                              dplyr::select(-gene) %>% unlist(),
                                            check.names=FALSE, stringsAsFactors = FALSE)
      }
    }
  }

  ## AR: Edge case when a gene only has values for single omic
  fix_index <- which(unlist(lapply(gene_tables_base, ncol)) == 1)
  if(length(fix_index)) {
    ## AR: account for having more than one gene with a single omic
    for(f in fix_index) {
      colnames(gene_tables_base[[f]]) <- paste0(names(gene_tables_base)[f], ".",
                                                colnames(gene_tables_base[[f]]))
    }
  }
  fix_index <- which(unlist(lapply(gene_tables_supp, ncol)) == 1)
  if(length(fix_index)) {
    ## AR: account for having more than one gene with a single omic
    for(f in fix_index) {
      colnames(gene_tables_supp[[f]]) <- paste0(names(gene_tables_supp)[f], ".",
                                                colnames(gene_tables_supp[[f]]))
    }

  }

  gene_tables_supp <- do.call("cbind", gene_tables_supp)
  gene_tables_base <- do.call("cbind", gene_tables_base)
  if(ncol(gene_tables_base) == nrow(base_data$rnaseq)) {
    colnames(gene_tables_base) <- base_data$rnaseq$gene
    if(ncol(gene_tables_supp) > 0)
      colnames(gene_tables_supp) <- base_data$rnaseq$gene
  }

  ## Remove elements with more missing elements than half the number of base individuals in either group
  remove_index <- unique(c(which(colSums(is.na(gene_tables_base)) > 0.5*nrow(base_data$rnaseq))#,
                           # which(colSums(is.na(gene_tables_supp)) > 0.5*nrow(supp_data$rnaseq)),
                           # which(apply(gene_tables_supp, 2, var) == 0),
                           # which(apply(gene_tables_base, 2, var) == 0)
  ))
  if(length(remove_index)) {
    gene_tables_supp <- gene_tables_supp[,-remove_index]
    gene_tables_base <- gene_tables_base[,-remove_index]
  }

  ## Only use missMDA for data imputation if specified by the user !!!
  ## Impute missing data for base and supp independently with missMDA (use type = "s" to scale data)
  lgr_tab <- strsplit(colnames(gene_tables_base), split = ".", fixed=TRUE) %>%
    lapply(function(x) {x[1]}) %>% unlist() %>% table()
  lgr <- as.numeric(lgr_tab)
  names(lgr) <- names(lgr_tab)
  if(sum(is.na(gene_tables_base)) & impute_MFA_missMDA) {
    gene_tables_base <- imputeMFA(gene_tables_base, group = lgr, ncp=2,
                                  type = rep("s", length(names(gene_tables_base))))$completeObs
  }
  if(sum(is.na(gene_tables_supp)) & impute_MFA_missMDA) {
    gene_tables_supp <- imputeMFA(gene_tables_supp, group = lgr, ncp=2,
                                  type = rep("s", length(names(gene_tables_supp))))$completeObs
  }

  gene_tables_base_s <- scale(gene_tables_base, center=TRUE, scale=TRUE)
  gene_tables_supp_s <-  t((t(gene_tables_supp)-attributes(gene_tables_base_s)$`scaled:center`) /
                             attributes(gene_tables_base_s)$`scaled:scale`)

  ## Combine all observations
  if(nrow(gene_tables_supp_s) > 0) {
    gene_tables <- rbind(gene_tables_base_s, gene_tables_supp_s)
  } else gene_tables <- gene_tables_base_s

  ## Run MFA: c = pre-scaled variables
  group <- c(rep("Base", nrow(base_data$clinical)), rep("Supp", nrow(supp_data$clinical)))
  if(nrow(gene_tables_supp_s) == 0) {
    ind.sup <- NULL
  } else {
    ind.sup=(nrow(gene_tables_base)+1):
      (nrow(gene_tables_base) + nrow(gene_tables_supp))
  }

  total_MFA <- MFA(gene_tables,
                   group=as.vector(lgr),
                   type=as.vector(rep("c",length(lgr))),
                   ind.sup=ind.sup,
                   ncp=ncol(gene_tables),
                   graph=FALSE,
                   name.group = names(lgr))

  ## Now calculate pathway deregulation score
  ps <- pathway_scores(total_MFA, ngenes=length(lgr))
  ps_df <- data.frame(bcr_patient_barcode = names(ps), group = group,
                      pathway_deregulation = ps, row.names=NULL, stringsAsFactors = FALSE)

  ## Now calculate the gene specific scores
  f <- total_MFA$ind$coord
  fk <- total_MFA$ind$coord.partiel
  gene_names <- gene_tables %>% colnames %>%
    strsplit(., ".", fixed=TRUE) %>% lapply(., function(x) x[1]) %>%
    unlist() %>% unique()
  fk_list <- vector("list", length(gene_names))
  names(fk_list) <- gene_names
  for(g in gene_names) {
    fk_list[[g]] <- fk[grep(paste0(".", g, "$"), rownames(fk)),]
  }
  df_final <- matrix(NA, nrow = nrow(f), ncol = length(gene_names))
  rownames(df_final) <- rownames(f)
  colnames(df_final) <- gene_names
  for(g in gene_names) {
    #   df_final[,g] <- rowSums(f*(fk_list[[g]] - f)) / rowSums(f^2)
    df_final[,g] <- rowSums((f*(fk_list[[g]] - f))) / sqrt(rowSums(f^2)) ## RENORMALIZED BY DISTANCE SCORE
  }
  gs_df <- df_final

  ## Omics: % contribution to each axis of the MFA
  omics_contrib_MFA <- NULL
  omics_groups <- strsplit(rownames(total_MFA$quanti.var$contrib), split = ".", fixed = TRUE) %>%
    lapply(function(xx) xx[2]) %>% unlist()
  if(!is.na(omics_groups[1])) {
    if(length(grep("hsa-mir", omics_groups))) {
      omics_groups[grep("hsa-mir", omics_groups)] <- "mirna"
    }
    omics_contrib_MFA <- round(rowsum(total_MFA$quanti.var$contrib, group = omics_groups), 2)
  }

  if(!full_results) {
    if(!is.na(omics_groups[1])) {
      omics_contrib_MFA_summary =
        data.frame(PC_10 = round(rowSums(omics_contrib_MFA[,1:min(10, ncol(omics_contrib_MFA))])/
                                   min(10, ncol(omics_contrib_MFA)), 2),
                   PC_all = round(rowSums(omics_contrib_MFA[,1:ncol(omics_contrib_MFA)])/
                                    ncol(omics_contrib_MFA), 2))
    } else {
      omics_contrib_MFA_summary <- data.frame(PC_10 = 1, PC_all = 1)
      row.names(omics_contrib_MFA_summary) <- "single-omic"
    }
  }

  if(full_results) {
    res <- list(Pathway_Name=Pathway_Name,
                pathway_deregulation = ps_df,
                pathway_gene_deregulation = gs_df,
                eig = total_MFA$eig,
                partial_coordinates_MFA = rbind(total_MFA$ind$coord.partiel,total_MFA$ind.sup$coord.partiel),
                ## Base individuals: % contributions to each PC
                ind_contrib_MFA = round(total_MFA$ind$contrib, 2),
                ## Genes: % contributions to each axis of the MFA
                gene_contrib_MFA = round(total_MFA$group$contrib, 2),
                ## Genes: Lg coefficients to show pairwise correlations among genes
                gene_Lg_MFA = round(total_MFA$group$Lg, 4),
                ## Omics: % contribution to each axis of the MFA
                omics_contrib_MFA = omics_contrib_MFA,
                ngenes = length(lgr),
                imputed_genes = imputed_genes,
                removed_genes = removed_genes,
                total_MFA = total_MFA,
                gene_tables = gene_tables)
  } else {
    res <- list(Pathway_Name=Pathway_Name,
                pathway_deregulation = ps_df,
                pathway_gene_deregulation = gs_df,
                eig = total_MFA$eig,
                ## Base individuals: % contributions to each first 10 PCs and all PCs
                ind_contrib_MFA_summary =
                  data.frame(PC_10 = round(rowSums(total_MFA$ind$contrib[,1:min(10, ncol(total_MFA$ind$contrib))])/
                                             min(10, ncol(total_MFA$ind$contrib)), 2),
                             PC_all = round(rowSums(total_MFA$ind$contrib[,1:ncol(total_MFA$ind$contrib)])/
                                              ncol(total_MFA$ind$contrib), 2)),
                ## Genes: % contributions to each first 10 PCs and all PCs
                gene_contrib_MFA_summary =
                  data.frame(PC_10 = round(rowSums(total_MFA$group$contrib[,1:min(10, ncol(total_MFA$group$contrib))])/
                                             min(10, ncol(total_MFA$group$contrib)), 2),
                             PC_all = round(rowSums(total_MFA$group$contrib[,1:ncol(total_MFA$group$contrib)])/
                                              ncol(total_MFA$group$contrib), 2)),
                ## Omics: % contribution to each first 10 PCs and all PCs
                omics_contrib_MFA_summary = omics_contrib_MFA_summary,
                ngenes = length(lgr),
                imputed_genes = imputed_genes,
                removed_genes = removed_genes)
  }
  class(res) <- "padma"
  return(res)
}



#-----------------------------------------------------------------------

## DO NOT EXPORT
pathway_scores <- function(MFA_results, ngenes) {
  return(deregulation=sqrt(rowSums(rbind(MFA_results$ind$coord,
                                         MFA_results$ind.sup$coord)^2)))/ngenes
}

## DO NOT EXPORT
pathway_subset <- function(x, y, subset, miRTarBase) {
  ## Initialize
  subset_rna_x <- subset_rna_y <- subset_methyl_x <- subset_methyl_y <-
    subset_cna_x <- subset_cna_y <- subset_mirna_x <- subset_mirna_y <- NULL

  if("rnaseq" %in% names(x) & "rnaseq" %in% names(y)) {
    subset_rna_x <- dplyr::filter(x$rnaseq, gene %in% subset)
    subset_rna_y <- dplyr::filter(y$rnaseq, gene %in% subset)
  }
  if("methyl" %in% names(x) & "methyl" %in% names(y)) {
    subset_methyl_x <- dplyr::filter(x$methyl, gene %in% subset)
    subset_methyl_y <- dplyr::filter(y$methyl, gene %in% subset)
  }
  if("mirna" %in% names(x) & "mirna" %in% names(y)) {
    miRTarBase_choose <- miRTarBase %>%
      filter(`Target Gene` %in% subset) %>%
      mutate(miRNA_lc = tolower(miRNA)) %>%
      dplyr::select(miRNA_lc, `Target Gene`) %>%
      filter(miRNA_lc %in% x$mirna$miRNA_lc) %>%
      unique()
    if(nrow(miRTarBase_choose) > 0) {
      subset_mirna_x <- data.frame(x$mirna, stringsAsFactors = FALSE, check.names=FALSE) %>%
        left_join(miRTarBase_choose, ., by = c("miRNA_lc", "Target Gene"))
      subset_mirna_y <- data.frame(y$mirna, stringsAsFactors = FALSE, check.names=FALSE) %>%
        left_join(miRTarBase_choose, ., by = c("miRNA_lc", "Target Gene"))
    } else { ## AR add: edge case where no miRs map to pathway genes
      subset_mirna_x  <- subset_mirna_y <- NULL
    }

  }
  if("cna" %in% names(x) & "cna" %in% names(y)) {
    subset_cna_x <- dplyr::filter(x$cna, gene %in% subset)
    subset_cna_y <- dplyr::filter(y$cna, gene %in% subset)
  }

  tmp <- list(x = list(clinical=x$clinical,
                       rnaseq=subset_rna_x,
                       methyl = subset_methyl_x,
                       mirna = subset_mirna_x,
                       cna = subset_cna_x),
              y = list(clinical=y$clinical,
                       rnaseq=subset_rna_y,
                       methyl = subset_methyl_y,
                       mirna = subset_mirna_y,
                       cna = subset_cna_y))

  return(tmp)
}


## DO NOT EXPORT
omics_subset <- function(x, index, apply_log) {
  xx <- vector("list", length = length(x))
  names(xx) <- names(x)
  for(i in names(x)) {
    ## Add one here because there is a gene column
    if(i == "clinical") {
      xx[[i]] <- x[[i]][index,,drop=FALSE]
    } else if(i != "mirna") {
      ## Apply log transformation as needed
      if(!is.null(apply_log)) {
        if(i %in% apply_log) {
          xx[[i]] <- cbind(x[[i]][,1, drop=FALSE], log(x[[i]][,c(index+1),drop=FALSE] + 1))
        } else {
          xx[[i]] <- x[[i]][,c(1,index+1), drop=FALSE]
        }
      } else {
        xx[[i]] <- x[[i]][,c(1,index+1), drop=FALSE]
      }
    } else {
      ## Apply log transformation as needed
      if(!is.null(apply_log)) {
        if(i %in% apply_log) {
          xx[[i]] <- cbind(x[[i]][,1:2, drop=FALSE], log(x[[i]][,c(index+2),drop=FALSE] + 1))
        } else {
          xx[[i]] <- x[[i]][,c(1:2,index+2), drop=FALSE]
        }
      } else {
        xx[[i]] <- x[[i]][,c(1:2,index+2), drop=FALSE]
      }
    }
  }
  return(xx)
}

## DO NOT EXPORT
whichminpositivevar <- function(xx) {
  xxx <- apply(xx, 1, var, na.rm=TRUE)
  return(which(xxx == min(xxx[xxx > 10e-5], na.rm=TRUE))[1])
}


