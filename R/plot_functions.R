
#' Plot an MFA factor map for individuals based on padma analysis
#'
#' Produce an MFA factor map for individuals for a pair of dimensions provided by the user, where
#' a normal confidence ellipsis is optionally superimposed.
#'
#' @param padma_obj Output from running the \code{padma} function
#' (with \code{"full_results = TRUE"})
#' @param dim_x Dimension number of the MFA to be plotted on the x-axis.
#' @param dim_y Dimension number of the MFA to be plotted on the y-axis.
#' @param plot_ellipse If \code{TRUE}, superimpose a normal confidence ellipsis
#' on the factor map.
#'
#' @return Factor map of class \code{ggplot}.
#' @export
#' @import ggplot2 ggrepel
plot_factor_map <- function(padma_obj, dim_x = 1, dim_y = 2, plot_ellipse = TRUE) {
  if(class(padma_obj) != "padma") stop("This plot function expects an object of class padma.")
  test <- padma_obj
  if(is.null(ncol(test$total_MFA$ind$coord)))
    stop("This plotting function may only be used if padma is run with full_results = TRUE.")
  if(dim_x >= ncol(test$total_MFA$ind$coord))
    stop("dim_x must be less than the largest dimension of the MFA.")
  if(dim_y >= ncol(test$total_MFA$ind$coord))
    stop("dim_y must be less than the largest dimension of the MFA.")


  axis_begin_1  <- floor(min(test$total_MFA$ind$coord[,dim_x]))
  axis_end_1    <- ceiling(max(test$total_MFA$ind$coord[,dim_x]))
  axis_begin_2  <- floor(min(test$total_MFA$ind$coord[,dim_y]))
  axis_end_2    <- ceiling(max(test$total_MFA$ind$coord[,dim_y]))
  lab_frame_1 <- data.frame(lab_1 = seq(axis_begin_1, axis_end_1, 2),
                            zero = 0) %>% subset(lab_1 != 0)
  lab_frame_2 <- data.frame(lab_2 = seq(axis_begin_2, axis_end_2, 2),
                            zero = 0) %>% subset(lab_2 != 0)
  # chart junk data
  tick_frame_1 <-
    data.frame(ticks_1 = seq(axis_begin_1, axis_end_1, by = 1),
               zero=0) %>% subset(ticks_1 != 0)
  tick_frame_2 <-
    data.frame(ticks_2 = seq(axis_begin_2, axis_end_2, by = 1),
               zero=0) %>% subset(ticks_2 != 0)
  tick_sz <- 0.1
  tmp <- data.frame(bcr_patient_barcode = rownames(test$total_MFA$ind$coord),
                    test$total_MFA$ind$coord[, c(dim_x, dim_y)],
                    stringsAsFactors = FALSE)
  colnames(tmp) <- c("bcr_patient_barcode", "Dim.1", "Dim.2")

  tmp$bcr_patient_barcode[-which(tmp$Dim.1 < quantile(tmp$Dim.1, prob = .025) |
                                   tmp$Dim.1 > quantile(tmp$Dim.1, prob = .975) |
                                   tmp$Dim.2 < quantile(tmp$Dim.2, prob = .025) |
                                   tmp$Dim.2 > quantile(tmp$Dim.2, prob = .975))] <- ""
  fig2b <- ggplot(tmp, aes(x=Dim.1,y=Dim.2)) +
    geom_segment(x = 0, xend = 0,
                 y = lab_frame_2$lab_2[1], yend = tail(lab_frame_2$lab_2, 1)+1,
                 size = 0.5) +
    geom_segment(y = 0, yend = 0,
                 x = lab_frame_1$lab_1[1], xend = tail(lab_frame_1$lab_1, 1)+1,
                 size = 0.5) +
    geom_segment(data = tick_frame_1,
                 aes(x = ticks_1, xend = ticks_1,
                     y = zero - tick_sz, yend = zero + tick_sz)) +
    geom_segment(data = tick_frame_2,
                 aes(x = zero - tick_sz, xend = zero + tick_sz,
                     y = ticks_2, yend = ticks_2))
  if(plot_ellipse == TRUE) {
    fig2b <- fig2b + stat_ellipse(color = 'grey90', alpha = 0.2, geom = "polygon",
                                  lwd = 1.25, fill ="grey90")
  }
  fig2b <- fig2b +
    geom_point(size=2, alpha = 0.6) +
    theme_minimal() +
    geom_text_repel(aes(label = bcr_patient_barcode)) +
    xlab(paste("MFA Dimension", dim_x)) + ylab(paste("MFA Dimension", dim_y))
  return(fig2b)
}





#' Plot an MFA partial factor map for a given individual and genes based on padma analysis
#'
#' Produce an MFA partial factor map for a given individual for a pair of dimensions provided
#' by the user.
#'
#' @param padma_obj Output from running the \code{padma} function
#' (with \code{"full_results = TRUE"})
#' @param id Index or sample name to be plotted. By default, the first sample is plotted.
#' @param dim_x Dimension number of the MFA to be plotted on the x-axis.
#' @param dim_y Dimension number of the MFA to be plotted on the y-axis.
#'
#' @return Partial factor map of class \code{ggplot}.
#' @export
#' @import ggplot2 ggrepel
plot_partial_factor_map <- function(padma_obj,
                                    id = padma_obj$pathway_deregulation$bcr_patient_barcode[1],
                                    dim_x = 1, dim_y = 2) {
  if(class(padma_obj) != "padma") stop("This plot function expects an object of class padma.")
  test <- padma_obj
  if(is.null(ncol(test$total_MFA$ind$coord)))
    stop("This plotting function may only be used if padma is run with full_results = TRUE.")
  if(dim_x >= ncol(test$total_MFA$ind$coord))
    stop("dim_x must be less than the largest dimension of the MFA.")
  if(dim_y >= ncol(test$total_MFA$ind$coord))
    stop("dim_y must be less than the largest dimension of the MFA.")

  ## Figure 3aa
  i <- id
  df <- data.frame(bcr_patient_barcode = rownames(test$total_MFA$ind$coord),
                   test$total_MFA$ind$coord) %>%
    gather(key = dimension, value = coord, -bcr_patient_barcode)
  pdf <- data.frame(bcr_patient_barcode = substr(rownames(test$total_MFA$ind$coord.partiel), 1, 12),
                    gene = unlist(lapply(strsplit(rownames(test$total_MFA$ind$coord.partiel), split = ".", fixed=TRUE), function(x) x[2])),
                    test$total_MFA$ind$coord.partiel) %>%
    gather(key = dimension, value = partiel_coord, -bcr_patient_barcode, -gene)
  df_choose <- df %>% filter(bcr_patient_barcode == i,
                             dimension %in% c(paste0("Dim.", dim_x),
                                              paste0("Dim.", dim_y))) %>%
    spread(key=dimension, value = coord)
  colnames(df_choose) <- c("bcr_patient_barcode", "Dim.1", "Dim.2")
  pdf_choose <- pdf %>% filter(bcr_patient_barcode == i,
                               dimension %in% c(paste0("Dim.", dim_x),
                                                paste0("Dim.", dim_y))) %>%
    spread(key=dimension, value = partiel_coord)
  colnames(pdf_choose) <- c("bcr_patient_barcode", "gene", "Dim.1", "Dim.2")
  pdf_choose$distance <- sqrt(rowSums(pdf_choose[,-c(1:2)]^2))
  axis_begin_1  <- floor(min(pdf_choose$Dim.1))
  axis_end_1    <- ceiling(max(pdf_choose$Dim.1))+1
  axis_begin_2  <- floor(min(pdf_choose$Dim.2))
  axis_end_2    <- ceiling(max(pdf_choose$Dim.2))+1
  lab_frame_1 <- data.frame(lab_1 = seq(axis_begin_1, axis_end_1, 2),
                            zero = 0) %>% subset(lab_1 != 0)
  lab_frame_2 <- data.frame(lab_2 = seq(axis_begin_2, axis_end_2, 2),
                            zero = 0) %>% subset(lab_2 != 0)
  tick_frame_1 <-
    data.frame(ticks_1 = seq(axis_begin_1, axis_end_1, by = 1),
               zero=0) %>% subset(ticks_1 != 0)
  tick_frame_2 <-
    data.frame(ticks_2 = seq(axis_begin_2, axis_end_2, by = 1),
               zero=0) %>% subset(ticks_2 != 0)
  tick_sz <- 0.1
  fig3aa <- ggplot(df_choose) +
    geom_point(aes(Dim.1, Dim.2), size =  5) +
    geom_point(data = pdf_choose, aes(Dim.1, Dim.2), alpha = 0.5) +
    geom_text_repel(data = pdf_choose, aes(Dim.1, Dim.2, label=gene)) +
    geom_segment(data = left_join(pdf_choose, df_choose, by ="bcr_patient_barcode"),
                 aes(x=Dim.1.x, xend = Dim.1.y, y = Dim.2.x, yend = Dim.2.y), alpha = 0.2, lty=2) +
    geom_segment(x = 0, xend = 0,
                 y = lab_frame_2$lab_2[1], yend = tail(lab_frame_2$lab_2, 1)+1,
                 size = 0.5) +
    geom_segment(y = 0, yend = 0,
                 x = lab_frame_1$lab_1[1], xend = tail(lab_frame_1$lab_1, 1)+1,
                 size = 0.5) +
    geom_segment(data = tick_frame_1,
                 aes(x = ticks_1, xend = ticks_1,
                     y = zero - tick_sz, yend = zero + tick_sz)) +
    geom_segment(data = tick_frame_2,
                 aes(x = zero - tick_sz, xend = zero + tick_sz,
                     y = ticks_2, yend = ticks_2)) +
    theme_minimal() +
    xlab(paste("MFA Dimension", dim_x)) + ylab(paste("MFA Dimension", dim_y))

  return(fig3aa)
}



#' Plot the omics contribution per MFA axis and the overall weighted contribution
#'
#' Plot barplots indicating the percent contribution of each omics to each MFA dimension,
#' as well as the overall weighted (by eigenvalue) percent contribution to the full
#' analysis.
#'
#' @param padma_obj Output from running the \code{padma} function
#' (with \code{"full_results = TRUE"})
#' @param max_dim Maxomum imension number of the MFA to be plotted.
#'
#' @return Barplots of percent variance contribution of class \code{ggplot}.
#' @export
#' @import ggplot2 ggrepel viridis cowplot
plot_omics_contrib <- function(padma_obj,
                                    max_dim = min(10, nrow(padma_obj$eig))) {
  if(class(padma_obj) != "padma") stop("This plot function expects an object of class padma.")
  test <- padma_obj
  if(is.null(ncol(test$total_MFA$ind$coord)))
    stop("This plotting function may only be used if padma is run with full_results = TRUE.")
  if(max_dim > ncol(test$total_MFA$ind$coord))
    stop("max_dim must be less than the largest dimension of the MFA.")

  ## Fig 5b: Omics contribution
  df <- data.frame(test$omics_contrib_MFA)
  df$Dim.0 <- apply(df, 1, weighted.mean,
                    w = test$eig[,"eigenvalue"] / sum(test$eig[,"eigenvalue"]))
  df$omics <- rownames(df)
  df <- df %>%
    gather(value = percent, key = dimension, -omics) %>%
    separate(dimension, into = c("var", "dimension")) %>%
    dplyr::select(-var) %>%
    mutate(dimension = as.numeric(dimension))
  df$dimension <- ifelse(df$dimension == 0, -5, df$dimension)
  eig_df <- data.frame(dimension = as.numeric(substr(row.names(test$eig), 6, 10)),
                       percentage_of_variance = as.numeric(test$eig[,2]))
  df <- full_join(df, eig_df, by = "dimension")
  g1 <- ggplot(filter(df, dimension <= max_dim, dimension > 0)) +
    geom_bar(aes(x = dimension, y = percent, fill = omics,
                 alpha = percentage_of_variance / 10),
             stat = "identity") +
    geom_text(data = filter(df, dimension <= max_dim, dimension > 0) %>% dplyr::select(dimension, percentage_of_variance) %>% unique(),
              aes(x = dimension, y = 105, label = round(percentage_of_variance, 2))) +
    scale_fill_viridis(discrete = TRUE) + scale_alpha(guide = 'none') +
    ylim(c(0, 105)) +
    theme_minimal_hgrid() +
    scale_x_continuous(breaks = 1:max_dim)
  g2 <- ggplot(filter(df, dimension == -5)) +
    geom_bar(aes(x = dimension, y = percent, fill = omics),
             stat = "identity") +
    scale_fill_viridis(discrete = TRUE, guide = "none") +
    xlab("Weighted average\nacross dimensions") +
    ylim(c(0, 105)) +
    theme_minimal_hgrid() +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank())
  fig5b <- plot_grid(g2, g1, nrow = 1, rel_widths = c(1,5), labels = c("", ""))
  return(fig5b)
}



