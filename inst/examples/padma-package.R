
LUAD_subset <- padma::LUAD_subset
## Create MultiAssayExperiment object with LUAD data
omics_data <- 
  list(rnaseq = LUAD_subset$rnaseq,
       methyl = LUAD_subset$methyl,
       mirna = LUAD_subset$mirna,
       cna = LUAD_subset$cna)
pheno_data <- 
  data.frame(LUAD_subset$clinical, 
             row.names = LUAD_subset$clinical$bcr_patient_barcode)
mae <-
  suppressMessages(
    MultiAssayExperiment::MultiAssayExperiment(
      experiments = omics_data, colData = pheno_data))

## Run padma
run_padma <- 
  padma(mae, gene_map = padma::mirtarbase,
        pathway_name = "c2_cp_BIOCARTA_D4GDI_PATHWAY", verbose = FALSE)

summary(run_padma)

## padma plots
factorMap(run_padma, dim_x = 1, dim_y = 2)
factorMap(run_padma, dim_x = 1, dim_y = 2,
           partial_id = "TCGA-78-7536")
omicsContrib(run_padma, max_dim = 10)
