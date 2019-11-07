
run_padma <- padma(LUAD_subset,
                   pathway_name = "c2_cp_BIOCARTA_D4GDI_PATHWAY")

plot_factor_map(run_padma, dim_x = 1, dim_y = 2)
plot_partial_factor_map(run_padma, dim_x = 1, dim_y = 2)
plot_partial_factor_map(run_padma, id = "TCGA-78-7536", dim_x = 1, dim_y = 2)
plot_omics_contrib(run_padma, max_dim = 10)
