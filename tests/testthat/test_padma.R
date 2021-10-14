LUAD_subset <- padma::LUAD_subset
omics_data <- 
  list(rnaseq = as.matrix(LUAD_subset$rnaseq),
       methyl = as.matrix(LUAD_subset$methyl),
       mirna = as.matrix(LUAD_subset$mirna),
       cna = as.matrix(LUAD_subset$cna))
pheno_data <- 
  data.frame(LUAD_subset$clinical, 
             row.names = LUAD_subset$clinical$bcr_patient_barcode)
mae <-
  suppressMessages(
    MultiAssayExperiment::MultiAssayExperiment(
      experiments = omics_data, colData = pheno_data))
run <- padma(mae, pathway_name = "c2_cp_BIOCARTA_D4GDI_PATHWAY",
             verbose = FALSE)
run_concise <- padma(mae, pathway_name = "c2_cp_BIOCARTA_D4GDI_PATHWAY",
             verbose = FALSE, full_results = FALSE)

##------------------------------------------------------------
context("Test of the function 'padma'.... testing inputs")

test_that("Function throws error if unexpected input", {
  expect_error(padma(as.matrix(do.call("cbind", LUAD_subset[-1])), 
                     pathway_name = "c2_cp_BIOCARTA_D4GDI_PATHWAY"))
  expect_error(padma(do.call("cbind", LUAD_subset[-1])), 
                     pathway_name = "c2_cp_BIOCARTA_D4GDI_PATHWAY")
})

test_that("Function throws error if phenotypes in list", {
  expect_error(padma(LUAD_subset, 
                     pathway_name = "c2_cp_BIOCARTA_D4GDI_PATHWAY"))
})

test_that("Function throws error if overlap in base_ids and supp_ids", {
  expect_error(padma(mae, 
                     base_ids = sampleMap(mae)$primary[1:10],
                     supp_ids = sampleMap(mae)$primary[5:15],
                     pathway_name = "c2_cp_BIOCARTA_D4GDI_PATHWAY"))
})

test_that("Function throws error if incorrect named pathway provided", {
  expect_error(padma(mae, 
                     pathway_name = "random_pathway"))
})

##------------------------------------------------------------
context("Test of the function 'padma'.... testing outputs")

test_that("Function returns empty results if no genes found", {
  expect_message(padma(mae, 
                     pathway_name = c("gene1", "gene2")))
  expect_null(padma(mae, pathway_name = c("gene1", "gene2")))
})

test_that("Function returns padmaResults object", {
  expect_s4_class(run, "padmaResults")
  expect_named(MFA_results(run), 
               c("eig","ind_contrib_MFA","gene_contrib_MFA",
                 "gene_Lg_MFA","omics_contrib_MFA", "total_MFA", 
                 "gene_tables"))
  expect_named(MFA_results(run_concise), 
               c("eig", "ind_contrib_MFA_summary", 
                 "gene_contrib_MFA_summary", "omics_contrib_MFA_summary"))
})

