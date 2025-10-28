#' Boxplots for selected features across groups
#'
#' @param X numeric matrix [features x samples].
#' @param group factor/character of length ncol(X) giving sample groups.
#' @param features character vector of feature (gene/protein) names to plot.
#' @param feature_order optional character vector to set facet order (defaults to `features` order).
#' @param palette optional named vector of group colors; if NULL, generated automatically.
#' @param ncol integer; facets per row (default 4).
#' @param add_jitter logical; add jittered points (default TRUE).
#' @param jitter_width numeric; jitter width (default 0.15).
#' @param notch logical; use notched boxplots (default FALSE).
#' @param free_y logical; facet with free y scales (default TRUE).
#' @param point_size numeric; jitter point size (default 1.6).
#' @param box_linewidth numeric; box outline width (default 0.4).
#' @param title,ylab,xlab plot labels.
#' @param save_path optional file path to save (pdf/png) using `ggplot2::ggsave`.
#' @param width,height numeric dimensions (inches) if saving.
#' @return ggplot object.
#' @export
plot_feature_boxplots <- function(
    X, group, features,
    feature_order = NULL,
    palette = NULL,
    ncol = 4,
    add_jitter = TRUE,
    jitter_width = 0.15,
    notch = FALSE,
    free_y = TRUE,
    point_size = 1.6,
    box_linewidth = 0.4,
    title = NULL,
    ylab = "log2 expression",
    xlab = NULL,
    save_path = NULL,
    width = 7, height = 4.5
) {
  stopifnot(is.matrix(X), length(group) == ncol(X))
  group <- factor(group)  # ensure factor for consistent legend order

  feats_present <- intersect(features, rownames(X))
  feats_missing <- setdiff(features, rownames(X))
  if (length(feats_present) == 0L) stop("None of the requested features are present in X.")
  if (length(feats_missing) > 0) {
    message("Skipping missing features: ", paste(feats_missing, collapse = ", "))
  }

  # long table
  mat <- X[feats_present, , drop = FALSE]
  df <- data.frame(
    Expression = as.numeric(mat),
    Gene = rep(rownames(mat), times = ncol(mat)),
    Group = rep(group, each = nrow(mat)),
    stringsAsFactors = FALSE
  )

  # facet order
  if (is.null(feature_order)) feature_order <- feats_present
  df$Gene <- factor(df$Gene, levels = feature_order)

  # auto palette if needed
  if (is.null(palette)) {
    levs <- levels(group)
    nlev <- length(levs)
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      pal <- RColorBrewer::brewer.pal(min(max(3, nlev), 8), "Set2")
      if (nlev > length(pal)) pal <- grDevices::rainbow(nlev)
    } else {
      pal <- grDevices::hcl.colors(nlev, "Set2", rev = FALSE)
    }
    palette <- setNames(pal[seq_len(nlev)], levs)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Group, y = Expression, fill = Group)) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, notch = notch, linewidth = box_linewidth, color = "black") +
    ggplot2::stat_boxplot(geom = "errorbar", width = 0.35, linewidth = box_linewidth, color = "black") +
    {if (add_jitter) ggplot2::geom_jitter(size = point_size, width = jitter_width, alpha = 0.6, shape = 21, color = "black", stroke = 0.2)} +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::facet_wrap(~ Gene, ncol = ncol, scales = if (free_y) "free_y" else "fixed") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::labs(title = title, y = ylab, x = xlab) +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      axis.text.x = ggplot2::element_text(size = 8, color = "black"),
      axis.text.y = ggplot2::element_text(size = 8, color = "black")
    )

  if (!is.null(save_path)) {
    ggplot2::ggsave(filename = save_path, plot = p, width = width, height = height, limitsize = FALSE)
  }
  p
}
