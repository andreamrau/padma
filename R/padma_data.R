#' MSigDB canonical pathways and corresponding gene lists
#'
#' List of 1322 pathways and corresponding gene symbols included in the
#' MSigDB canonical pathways curated gene set catalog, which
#' includes genes whose products are involved in metabolic and
#' signaling pathways reported in curated public databases. This
#' specifically corresponds to the "C2 curated gene sets" catalog
#' from MSigDB v5.2 available at
#' http://bioinf.wehi.edu.au/software/MSigDB/ as described in the
#' limma Bioconductor package.
#'
#' @docType data
#'
#' @usage data(msigdb)
#'
#' @format An object of class \code{data.frame} with two
#' columns: \code{geneset}, which provides the 1322 MSigDB curated
#' pathway names (e.g., \code{"c2_cp_BIOCARTA_41BB_PATHWAY"}) and
#' \code{symbol}, which provides the comma-separated corresponding
#' list of gene symbols.
#'
#' @keywords datasets
#'
#' @references
#' Liberzon et al. (2011) Bioinformatics 27:12, 1739-1740.
#' \href{https://doi.org/10.1093/bioinformatics/btr260}.
#'
#'
#' @source
#' \href{http://software.broadinstitute.org/gsea/msigdb/genesets.jsp?collection=CP}{MSigDB Gene sets}
#' \href{http://bioinf.wehi.edu.au/software/MSigDB/}
#'
#' @examples
#' data(msigdb)
#' head(msigdb)
"msigdb"


#' Curated miR-target interaction predictions from miRTarBase
#'
#' List of 10,754 predicted miRNA gene targets from miRTarBase
#' (version 7.0), filtered to include only predictions with the
#' "Functional MTI" support type.
#'
#' @docType data
#'
#' @usage data(mirtarbase)
#'
#' @format An object of class \code{data.frame} with two
#' columns: \code{miRNA}, which provides the miRNA identifier
#' (e.g., \code{"hsa-miR-20a-5p"}) and
#' \code{Target Gene}, which provides the corresponding
#' predicted gene target.
#'
#' @keywords datasets
#'
#' @references
#' Chou et al. (2018) Nucleic Acids Research 46, D296–D302.
#' \href{https://doi.org/10.1093/nar/gkx1067}.
#'
#' @source
#' \href{http://mirtarbase.mbc.nctu.edu.tw/php/index.php}
#'
#' @examples
#' data(mirtarbase)
#' head(mirtarbase)
"mirtarbase"


#' Subset of batch-corrected multi-omic TCGA data in lung adenocarcinoma
#'
#' Multi-omic (RNA-seq, copy number alterations, methylation, and
#' miRNA-seq) and phenotypic data in 144 individuals in the TCGA-LUAD
#' data for the 13 genes in the D4-GDI signaling pathway.
#'
#' @docType data
#'
#' @usage data(LUAD_subset)
#'
#' @format A named list of five objects of class \code{data.frame}
#' containing a subset of the batch-corrected multi-omic TCGA
#' data from lung adenocarcinoma, corresponding to the 13 genes
#' in the D4 GDI signaling pathway:
#' \code{"clinical"} is of dimension 144 x 55 and contains clinical
#' variables for the 144 individuals. \code{"rnaseq"}, \code{"methyl"},
#' and \code{"cna"} are of dimension 13 (genes) x 145, where the
#' first column \code{"gene"} contains the gene symbol and the
#' remaining 144 columns correspond to the TCGA individuals.
#' \code{"mirna"} is of dimension 1 x 146 for the single miRNA
#' predicted to target one of the D4 GDI signaling pathway genes,
#' and the first two columns \code{miRNA} and \code{Target Gene}
#' respectively provide the miRNA identifier
#' (\code{"hsa-mir-421"}) and the corresponding
#' predicted gene target (\code{"hsa-mir-421"}).
#'
#' @keywords datasets
#'
#' @references
#' The Cancer Genome Atlas Research Network (2014) Nature 511, 543-550.
#' \href{https://doi.org/10.1038/nature13385}.
#'
#' Rau et al. (2019) bioRxiv, \href{https://doi.org/10.1101/827022}.
#'
#'
#' @source
#' \href{https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga}{TCGA}
#'
#' @examples
#' data(LUAD_subset)
#' head(LUAD_subset)
"LUAD_subset"

