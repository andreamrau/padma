## ----setup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- fig.align = "center", out.width = "50%", echo=FALSE------------------
knitr::include_graphics("hex_padma_v2.png")

## --------------------------------------------------------------------------
library(padma)

## ---- eval = FALSE---------------------------------------------------------
#  D4GDI <- msigdb[grep("D4GDI", msigdb$geneset),
#                  "geneset"]
#  run_padma <- padma(LUAD_subset,
#                     pathway_name = D4GDI)

## ---- eval = FALSE---------------------------------------------------------
#  plot_factor_map(run_padma)
#  plot_partial_factor_map(run_padma,
#                          id = "TCGA-78-7536")
#  plot_omics_contrib(run_padma)

## ---- explore_LUAD_subset--------------------------------------------------
names(LUAD_subset)

lapply(LUAD_subset, class)

lapply(LUAD_subset, dim)

## ---- msigdb---------------------------------------------------------------
head(msigdb)

## --------------------------------------------------------------------------
head(mirtarbase)

## ---- runpadma-------------------------------------------------------------
D4GDI <- msigdb[grep("D4GDI", msigdb$geneset), "geneset"]
run_padma <- padma(LUAD_subset,
                   pathway_name = D4GDI)

## ---- runpadma2, eval=FALSE------------------------------------------------
#  D4GDI_genes <- unlist(strsplit(msigdb[grep("D4GDI", msigdb$geneset), "symbol"],
#                          ", "))
#  D4GDI_genes
#  run_padma_again <- padma(LUAD_subset,
#                     pathway_name = D4GDI_genes)

## --------------------------------------------------------------------------
plot_factor_map(run_padma, dim_x = 1, dim_y = 2)

## --------------------------------------------------------------------------
plot_partial_factor_map(run_padma,
                        id = "TCGA-78-7536",
                        dim_x = 1,
                        dim_y = 2)

## --------------------------------------------------------------------------
plot_omics_contrib(run_padma, max_dim = 10)

## ---- eval = FALSE---------------------------------------------------------
#  padma(LUAD_subset,
#        pathway_name = D4GDI,
#        base_ids = 1:10,
#        supp_ids = 15:20, ...length())

## ---- eval = FALSE---------------------------------------------------------
#  padma(LUAD_subset,
#        pathway_name = D4GDI,
#        apply_log = "rnaseq", ...)

## ---- eval = FALSE---------------------------------------------------------
#  padma(LUAD_subset,
#        pathway_name = D4GDI,
#        impute_MFA_missMDA = TRUE, ...)

## ---- eval = FALSE---------------------------------------------------------
#  padma(LUAD_subset,
#        pathway_name = D4GDI,
#        full_results = FALSE, ...)

