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
context("Test of 'padma' plotting.... testing factorMap inputs")

test_that("Function throws error if full_results = FALSE", {
  expect_error(factorMap(run_concise))
  expect_error(factorMap(run_concise, partial_id = "TCGA-78-7536"))
})

test_that("Function throws error if not a padmaResults objet", {
  expect_error(factorMap(mae))
})

test_that("Function throws error if dim_x or dim_y not appropriate values", {
  expect_error(factorMap(run, dim_x = -10, dim_y = 2))
  expect_error(factorMap(run, dim_x = 1000, dim_y = 2000))
})

test_that("Function throws error if ggrepel = TRUE but ggplot = FALSE", {
  expect_error(factorMap(run, ggrepel = TRUE, ggplot = FALSE))
})

test_that("Function throws error if incorrect sample id provided", {
  expect_error(factorMap(run, partial_id = "fake_sample"))
})

##------------------------------------------------------------
context("Test of 'padma' plotting.... testing factorMap outputs")

test_that("Function returns ggplot if ggplot = TRUE", {
  expect_is(factorMap(run, ggplot = TRUE), "ggplot")
})

##------------------------------------------------------------
context("Test of 'padma' plotting.... testing omicsContrib inputs")


test_that("Function throws error if not a padmaResults objet", {
  expect_error(omicsContrib(mae))
})

test_that("Function throws error if full_results = FALSE", {
  expect_error(omicsContrib(run_concise))
})

test_that("Function throws error if max_dim not appropriate values", {
  expect_error(omicsContrib(run, max_dim = -1000))
  expect_error(omicsContrib(run, max_dim = 1000))
})

##------------------------------------------------------------
context("Test of 'padma' plotting.... testing omicsContrib outputs")

test_that("Function returns ggplot if ggplot = TRUE", {
  expect_is(omicsContrib(run, ggplot = TRUE), "ggplot")
})

