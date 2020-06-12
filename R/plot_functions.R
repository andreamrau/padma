#' Plot an MFA factor map for individuals or partial factor map based on 
#' padma analysis
#'
#' Produce an MFA factor map for individuals, or MFA partial factor map for
#' a given individual, for a pair of dimensions provided by the user.
#'
#' @param padma_obj Output from running the \code{padma} function
#' (with \code{'full_results = TRUE'})
#' @param partial_id Index or sample name to be plotted for a 
#' partial factor map.
#' @param dim_x Dimension number of the MFA to be plotted on the x-axis.
#' @param dim_y Dimension number of the MFA to be plotted on the y-axis.
#' @param plot_ellipse If \code{TRUE}, superimpose a normal confidence ellipsis
#' on the factor map.
#' @param ggplot If \code{TRUE}, use \code{ggplot2} for plotting
#' @param repel_labels If \code{TRUE}, use \code{ggrepel} to repel sample 
#' labels from each other
#'
#' @return Plot, or factor map of class \code{ggplot} if \code{ggplot2 = TRUE}.
#' @export
#' @example inst/examples/padma-package.R
#' @importFrom graphics abline barplot layout lines mtext par points text
factorMap <- function(padma_obj, partial_id = NULL, dim_x = 1, dim_y = 2, 
                      plot_ellipse = TRUE, ggplot = TRUE, 
                      repel_labels = ifelse(ggplot == TRUE, TRUE, FALSE)) {
    if (!is(padma_obj, "padmaResults")) 
        stop("This plot function expects an object of class padmaResults.")
    test <- padma_obj
    if (is.null(MFA_results(test)$total_MFA)) 
 stop("Plotting may only be used if padma is run with full_results = TRUE.")
    if (dim_x >= ncol(MFA_results(test)$total_MFA$ind$coord)) 
        stop("dim_x must be less than the largest dimension of the MFA.")
    if (dim_y >= ncol(MFA_results(test)$total_MFA$ind$coord)) 
        stop("dim_y must be less than the largest dimension of the MFA.")
    if (repel_labels & !ggplot) 
        stop("ggplot plotting must be used if repel_labels = TRUE.")
    
    lab_1 <- lab_2 <- ticks_1 <- ticks_2 <- NULL
    
    ## Full factor map -------------------------------------------------
    if (is.null(partial_id)) {
        axis_begin_1 <- 
            floor(min(MFA_results(test)$total_MFA$ind$coord[, dim_x]))
        axis_end_1 <- 
            ceiling(max(MFA_results(test)$total_MFA$ind$coord[, dim_x]))
        axis_begin_2 <- 
            floor(min(MFA_results(test)$total_MFA$ind$coord[, dim_y]))
        axis_end_2 <- 
            ceiling(max(MFA_results(test)$total_MFA$ind$coord[, dim_y]))
        lab_frame_1 <- subset(data.frame(lab_1 = seq(axis_begin_1, axis_end_1, 
            2), zero = 0), lab_1 != 0)
        lab_frame_2 <- subset(data.frame(lab_2 = seq(axis_begin_2, axis_end_2, 
            2), zero = 0), lab_2 != 0)
        # chart junk data
        tick_frame_1 <- subset(data.frame(ticks_1 = seq(axis_begin_1, 
                                                        axis_end_1, by = 1), 
                                          zero = 0), ticks_1 != 0)
        tick_frame_2 <- subset(data.frame(ticks_2 = seq(axis_begin_2, 
                                                        axis_end_2, by = 1), 
                                          zero = 0), ticks_2 != 0)
        tick_sz <- 0.1
        tick_frame_1$m_tick_sz <- tick_frame_1$zero - tick_sz
        tick_frame_1$p_tick_sz <- tick_frame_1$zero + tick_sz
        tick_frame_2$m_tick_sz <- tick_frame_2$zero - tick_sz
        tick_frame_2$p_tick_sz <- tick_frame_2$zero + tick_sz
        
        tmp <- data.frame(
            primary = rownames(MFA_results(test)$total_MFA$ind$coord), 
            MFA_results(test)$total_MFA$ind$coord[, c(dim_x, dim_y)], 
            stringsAsFactors = FALSE)
        colnames(tmp) <- c("primary", "Dim.1", "Dim.2")
        
        ## Remove all but most extreme sample labels
        tmp_index <- which(tmp$Dim.1 < quantile(tmp$Dim.1, prob = 0.025) | 
                               tmp$Dim.1 > quantile(tmp$Dim.1, prob = 0.975) | 
                               tmp$Dim.2 < quantile(tmp$Dim.2, prob = 0.025) | 
                               tmp$Dim.2 > quantile(tmp$Dim.2, prob = 0.975))
        tmp$primary[-tmp_index] <- ""
        
        ## Plot with ggplot
        if (ggplot) {
            if (!requireNamespace("ggplot2")) 
                stop("Please load ggplot2 to plot with ggplot.")
            fig2b <- ggplot2::ggplot(tmp, 
                                     ggplot2::aes_string(x = "Dim.1", 
                                                         y = "Dim.2")) + 
                ggplot2::geom_segment(x = 0, xend = 0, 
                                      y = lab_frame_2$lab_2[1], 
                  yend = tail(lab_frame_2$lab_2, 1) + 1, size = 0.5) + 
                ggplot2::geom_segment(y = 0, yend = 0, 
                                      x = lab_frame_1$lab_1[1], 
                  xend = tail(lab_frame_1$lab_1, 1) + 1, size = 0.5) + 
                ggplot2::geom_segment(data = tick_frame_1, 
                  ggplot2::aes_string(x = "ticks_1", xend = "ticks_1", 
                                      y = "m_tick_sz", yend = "p_tick_sz")) + 
                ggplot2::geom_segment(data = tick_frame_2, 
                  ggplot2::aes_string(x = "m_tick_sz", xend = "p_tick_sz", 
                                      y = "ticks_2", 
                  yend = "ticks_2"))
            if (plot_ellipse) {
                fig2b <- fig2b + 
                    ggplot2::stat_ellipse(color = "grey90", alpha = 0.2, 
                                          geom = "polygon", lwd = 1.25, 
                                          fill = "grey90")
            }
            fig2b <- fig2b + 
                ggplot2::geom_point(size = 2, alpha = 0.6) + 
                ggplot2::theme_minimal() + 
                ggplot2::xlab(paste("MFA Dimension", dim_x)) + 
                ggplot2::ylab(paste("MFA Dimension", dim_y))
            if (repel_labels) {
                if (!requireNamespace("ggrepel")) 
                  stop("Please load ggrepel to repel text labels.")
                fig2b <- fig2b + 
                    ggrepel::geom_text_repel(
                        ggplot2::aes_string(label = "primary"))
            } else {
                fig2b <- fig2b + 
                    ggplot2::geom_text(ggplot2::aes_string(label = "primary"))
            }
            return(fig2b)
        } else {
            if (plot_ellipse) {
                if (!requireNamespace("car")) 
                  stop("Please load car to draw ellipse in base graphics.")
                car::dataEllipse(as.matrix(tmp[, 2:3]), level = 0.95, 
                                 robust = TRUE, 
                  center.pch = 0, xlab = paste("MFA Dimension", dim_x), 
                  ylab = paste("MFA Dimension", 
                    dim_y), xlim = c(axis_begin_1, axis_end_1), 
                  ylim = c(axis_begin_2, 
                    axis_end_2))
                text(x = tmp[, 2], y = tmp[, 3], label = tmp[, 1])
                abline(h = 0, lty = 2)
                abline(v = 0, lty = 2)
            } else {
                plot(x = tmp[, 2], y = tmp[, 3], xlab = paste("MFA Dimension", 
                  dim_x), ylab = paste("MFA Dimension", dim_y), 
                  xlim = c(axis_begin_1, 
                  axis_end_1), ylim = c(axis_begin_2, axis_end_2))
                text(x = tmp[, 2], y = tmp[, 3], label = tmp[, 1])
                abline(h = 0, lty = 2)
                abline(v = 0, lty = 2)
            }
        }
    } else {
        ## Partial factor map -------------------------------------------------
        i <- partial_id
        tmp <- MFA_results(test)$total_MFA$ind$coord
        df <- data.frame(rbind(cbind(rownames(tmp), colnames(tmp)[dim_x], tmp[, 
            dim_x]), cbind(rownames(tmp), colnames(tmp)[dim_y], tmp[, dim_y])), 
            check.rows = FALSE)
        colnames(df) <- c("primary", "dimension", "coord")
        
        tmp2 <- data.frame(
            primary = substr(rownames(
                MFA_results(test)$total_MFA$ind$coord.partiel), 
            1, 12), 
            gene = unlist(lapply(strsplit(
                rownames(MFA_results(test)$total_MFA$ind$coord.partiel), 
            split = ".", fixed = TRUE), 
            function(x) x[2])), 
            MFA_results(test)$total_MFA$ind$coord.partiel)
        index <- which(colnames(tmp2) %in% 
                           c(paste0("Dim.", dim_x), paste0("Dim.", dim_y)))
        pdf_choose <- tmp2[which(tmp2$primary == i), c(1, 2, index)]
        colnames(pdf_choose) <- c("primary", "gene", "Dim.1", "Dim.2")
        df_choose <- df[which(df$primary == i & df$dimension %in% 
                                  c(paste0("Dim.", dim_x), 
                                    paste0("Dim.", dim_y))), ]
        df_choose <- data.frame(primary = df_choose$primary[1], 
                                as.numeric(df_choose$coord[1]), 
                                as.numeric(df_choose$coord[2]))
        colnames(df_choose) <- c("primary", "Dim.1", "Dim.2")
        
        if(!nrow(pdf_choose))
          stop("Sample ID not found: please check the partial_id argument.")
        pdf_choose$distance <- sqrt(rowSums(pdf_choose[, -c(1, 2)]^2))
        axis_begin_1 <- floor(min(pdf_choose[, 3]))
        axis_end_1 <- ceiling(max(pdf_choose[, 3])) + 1
        axis_begin_2 <- floor(min(pdf_choose[, 4]))
        axis_end_2 <- ceiling(max(pdf_choose[, 4])) + 1
        lab_frame_1 <- subset(data.frame(lab_1 = seq(axis_begin_1, axis_end_1, 
            2), zero = 0), lab_1 != 0)
        lab_frame_2 <- subset(data.frame(lab_2 = seq(axis_begin_2, axis_end_2, 
            2), zero = 0), lab_2 != 0)
        tick_frame_1 <- subset(data.frame(ticks_1 = seq(axis_begin_1, 
                                                        axis_end_1, 
            by = 1), zero = 0), ticks_1 != 0)
        tick_frame_2 <- subset(data.frame(ticks_2 = seq(axis_begin_2, 
                                                        axis_end_2, 
            by = 1), zero = 0), ticks_2 != 0)
        tick_sz <- 0.1
        
        tick_frame_1$m_tick_sz <- tick_frame_1$zero - tick_sz
        tick_frame_1$p_tick_sz <- tick_frame_1$zero + tick_sz
        tick_frame_2$m_tick_sz <- tick_frame_2$zero - tick_sz
        tick_frame_2$p_tick_sz <- tick_frame_2$zero + tick_sz
        
        ## Plot with ggplot
        if (ggplot) {
            if (!requireNamespace("ggplot2")) 
                stop("Please load ggplot2 to plot with ggplot.")
            fig3aa <- ggplot2::ggplot(df_choose) + 
                ggplot2::geom_point(ggplot2::aes_string("Dim.1", 
                  "Dim.2"), size = 5) + 
                ggplot2::geom_point(data = pdf_choose, 
                  ggplot2::aes_string("Dim.1", 
                "Dim.2"), alpha = 0.5) + 
                ggplot2::geom_segment(
                  data = merge(pdf_choose,  df_choose, by = "primary"), 
                  ggplot2::aes_string(x = "Dim.1.x", xend = "Dim.1.y", 
                                      y = "Dim.2.x", yend = "Dim.2.y"), 
                  alpha = 0.2, lty = 2) + 
                ggplot2::geom_segment(x = 0, xend = 0, 
                                      y = lab_frame_2$lab_2[1], 
                  yend = tail(lab_frame_2$lab_2, 1) + 1, size = 0.5) + 
                ggplot2::geom_segment(y = 0, 
                  yend = 0, x = lab_frame_1$lab_1[1], 
                  xend = tail(lab_frame_1$lab_1, 1) + 1, size = 0.5) + 
                ggplot2::geom_segment(data = tick_frame_1, 
                  ggplot2::aes_string(x = "ticks_1", xend = "ticks_1", 
                                      y = "m_tick_sz", yend = "p_tick_sz")) + 
                ggplot2::geom_segment(data = tick_frame_2, 
                  ggplot2::aes_string(x = "m_tick_sz", xend = "p_tick_sz", 
                                      y = "ticks_2",  yend = "ticks_2")) + 
                ggplot2::theme_minimal() + 
                ggplot2::xlab(paste("MFA Dimension", dim_x)) + 
                ggplot2::ylab(paste("MFA Dimension", dim_y))
            if (repel_labels) {
                if (!requireNamespace("ggrepel")) 
                  stop("Please load ggrepel to repel text labels.")
                fig3aa <- fig3aa + ggrepel::geom_text_repel(data = pdf_choose, 
                  ggplot2::aes_string("Dim.1", "Dim.2", label = "gene"))
                
            } else {
                fig3aa <- fig3aa + 
                    ggplot2::geom_text(data = pdf_choose, 
                                       ggplot2::aes_string("Dim.1", "Dim.2", 
                                                           label = "gene"))
            }
            return(fig3aa)
        } else {
            ## Plot with base
            plot(df_choose$Dim.1, df_choose$Dim.2, 
                 xlab = paste("MFA Dimension", 
                dim_x), ylab = paste("MFA Dimension", dim_y), 
                xlim = c(axis_begin_1, 
                axis_end_1), ylim = c(axis_begin_2, axis_end_2))
            abline(h = 0, lty = 1)
            abline(v = 0, lty = 1)
            for (ii in seq_len(nrow(pdf_choose))) {
                lines(rbind(pdf_choose[ii, 3:4], df_choose[, -1]), lty = 2, 
                      col = "gray60")
                points(pdf_choose[ii, 3:4], lty = 2, col = "gray60", pch = 20)
                text(pdf_choose[ii, 3:4], label = pdf_choose$gene[ii])
            }
            points(df_choose$Dim.1, df_choose$Dim.2, cex = 3, pch = 20)
        }
    }
}





#' Plot the omics contribution per MFA axis and the overall weighted 
#' contribution
#'
#' Plot barplots indicating the percent contribution of each omics to each MFA 
#' dimension, as well as the overall weighted (by eigenvalue) percent 
#' contribution to the full analysis.
#'
#' @param padma_obj Output from running the \code{padma} function
#' (with \code{'full_results = TRUE'})
#' @param max_dim Maximum dimension number of the MFA to be plotted
#' @param ggplot If \code{TRUE}, use \code{ggplot2} for plotting (and
#' \code{cowplot} for combining ggplots)
#' @return Barplots of percent variance contribution, optionally of class 
#' \code{ggplot}.
#' @export
#' @example inst/examples/padma-package.R
#' @importFrom stats na.omit quantile reshape weighted.mean
#' @importFrom utils tail
omicsContrib <- function(padma_obj, 
                         max_dim = min(10, nrow(MFA_results(padma_obj)$eig)), 
                         ggplot = TRUE) {
    
    if (!is(padma_obj, "padmaResults")) 
        stop("This plot function expects an object of class padmaResults.")
    test <- padma_obj
    if (is.null(MFA_results(test)$total_MFA)) 
  stop("Plotting may only be used if padma is run with full_results = TRUE.")
    if (max_dim > ncol(MFA_results(test)$total_MFA$ind$coord)) 
        stop("max_dim must be less than the largest dimension of the MFA.")
    
    df <- data.frame(MFA_results(test)$omics_contrib_MFA)
    df$Dim.0 <- apply(df, 1, weighted.mean, 
                      w = MFA_results(test)$eig[, "eigenvalue"] / 
                          sum(MFA_results(test)$eig[, 
        "eigenvalue"]))
    df$omics <- rownames(df)
    
    df <- reshape(df, varying = -ncol(df), timevar = "dimension", 
                  v.names = "percent", 
        direction = "long", sep = ".", 
        times = as.numeric(unlist(lapply(strsplit(colnames(df)[-ncol(df)], 
            split = ".", fixed = TRUE), function(x) x[2]))))
    df$dimension <- ifelse(df$dimension == 0, -5, df$dimension)
    
    eig_df <- data.frame(
        dimension = as.numeric(substr(row.names(
          MFA_results(test)$eig), 6, 10)), 
        percentage_of_variance = as.numeric(MFA_results(test)$eig[, 2]))
    df <- merge(df, eig_df, by = "dimension", all.x = TRUE)
    df$percentage_of_variance_10 <- df$percentage_of_variance/10
    df$percentage_of_variance_round <- round(df$percentage_of_variance, 2)
    
    if (ggplot) {
        ## Plot with ggplot
        if (!requireNamespace("ggplot2") | !requireNamespace("cowplot")) 
            stop("Please load ggplot2 and cowplot to plot with ggplot.")
        
        g1 <- ggplot2::ggplot(df[which(df$dimension <= max_dim & 
                                           df$dimension > 0), ]) + 
            ggplot2::geom_bar(ggplot2::aes_string(x = "dimension", 
                                                  y = "percent", 
              fill = "omics", alpha = "percentage_of_variance_10"), 
              stat = "identity") + 
            ggplot2::geom_text(data = unique(df[which(df$dimension <= max_dim & 
                df$dimension > 0), c("dimension", "percentage_of_variance", 
                                     "percentage_of_variance_round")]), 
              ggplot2::aes_string(x = "dimension", y = 105, 
                                  label = "percentage_of_variance_round")) + 
            ggplot2::scale_alpha(guide = "none") + ggplot2::ylim(c(0, 105)) + 
            cowplot::theme_minimal_hgrid() + 
            ggplot2::scale_x_continuous(breaks = seq_len(max_dim)) + 
            ggplot2::ylab("")
        g2 <- ggplot2::ggplot(df[which(df$dimension == -5), ]) + 
            ggplot2::geom_bar(ggplot2::aes_string(x = "dimension", 
            y = "percent", fill = "omics"), stat = "identity") + 
            ggplot2::xlab("Weighted average\nacross dimensions") + 
            ggplot2::ylim(c(0, 105)) + cowplot::theme_minimal_hgrid() + 
            ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
            axis.text.x = ggplot2::element_blank(), legend.position = "none")
        fig5b <- cowplot::plot_grid(g2, g1, nrow = 1, rel_widths = c(1, 5), 
                                    labels = c("", 
            ""))
        return(fig5b)
    } else {
        ## Plot with base:
        op <- par()
        layout(matrix(c(1, 2, 2, 2), nrow = 1))
        barplot(as.matrix(df[which(df$dimension == -5), c("percent"), 
                             drop = FALSE]), 
            names.arg = "", xlab = "Weighted average\nacross dimensions", 
            ylab = "percent")
        df2 <- na.omit(df[which(df$dimension <= max_dim & df$dimension > 0), 
                          c("omics", 
            "percent", "dimension")])
        df2 <- reshape(df2, v.names = "percent", idvar = "omics", 
                       timevar = "dimension", 
            direction = "wide")
        rownames(df2) <- df2[, 1]
        df2 <- df2[, -1]
        b <- barplot(as.matrix(df2), 
                     names.arg = unlist(lapply(strsplit(colnames(df2), 
            split = ".", fixed = TRUE), function(x) x[2])), xlab = "Dimension", 
            legend.text = TRUE, args.legend = list(horiz = TRUE, bty = "n", 
                                                   x = ncol(df2), 
                y = 115))
        mtext(at = b, side = 3, round(unique(df[which(df$dimension <= max_dim & 
            df$dimension > 0), "percentage_of_variance"]), 1))
        suppressWarnings(par(op))
    }
}
