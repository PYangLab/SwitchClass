#' Plot heatmap for quadrant-classified features
#'
#' Generates a standardized heatmap of feature expression grouped by quadrant
#' (Q1–Q4) and sample group (e.g., HC / Pre / Post). Quadrants are optionally
#' filtered by delta cutoff thresholds.
#'
#' @param X numeric matrix [features x samples]; expression matrix.
#' @param group factor/character of length = ncol(X); sample group labels.
#' @param df_quad data.frame with columns: \code{feature}, \code{quadrant}, and \code{delta}.
#' @param cutoff numeric; delta threshold for selecting strong quadrant features (default 0.15).
#'   Quadrants Q1 and Q4 use \code{delta > cutoff}; Q2 and Q3 use \code{delta < -cutoff}.
#' @param palette_group optional named color vector for sample groups.
#' @param palette_quad optional named color vector for quadrants.
#' @param scale logical; z-score features (default TRUE).
#' @param cap numeric; cap scaled values at ±cap (default 2).
#' @param show_rownames,show_colnames logical; show labels (default TRUE).
#' @param cluster_rows,cluster_cols logical; cluster rows/columns (default FALSE).
#' @param width,height numeric; optional dimensions for PDF export.
#' @return Invisibly returns the \code{pheatmap} object.
#' @examples
#' \dontrun{
#' plot_heatmap_quadrants_auto(
#'   X = crc.imputed,
#'   group = group,
#'   df_quad = df,
#'   cutoff = 0.15
#' )
#' }
#' @export
plot_heatmap_quadrants <- function(
    X, group, df_quad,
    cutoff = 0.15,
    palette_group = NULL,
    palette_quad  = NULL,
    scale = TRUE, cap = 2,
    show_rownames = TRUE, show_colnames = TRUE,
    cluster_rows = FALSE, cluster_cols = FALSE
) {
  if (!requireNamespace("pheatmap", quietly = TRUE))
    stop("Package 'pheatmap' is required.")
  stopifnot(is.matrix(X), length(group) == ncol(X))
  stopifnot(all(c("feature","quadrant","delta") %in% names(df_quad)))

  # --- Select features per quadrant based on cutoff ---
  g1 <- df_quad$feature[df_quad$quadrant == "Q1" & df_quad$delta >  cutoff]
  g2 <- df_quad$feature[df_quad$quadrant == "Q2" & df_quad$delta < -cutoff]
  g3 <- df_quad$feature[df_quad$quadrant == "Q3" & df_quad$delta < -cutoff]
  g4 <- df_quad$feature[df_quad$quadrant == "Q4" & df_quad$delta >  cutoff]
  quad_genes <- list(Q1 = g1, Q2 = g2, Q3 = g3, Q4 = g4)

  # --- Combine valid features ---
  genes_ordered <- unlist(quad_genes)
  genes_ordered <- intersect(genes_ordered, rownames(X))
  if (length(genes_ordered) == 0)
    stop("No quadrant features found above cutoff.")
  mat <- X[genes_ordered, , drop = FALSE]

  # --- Z-score and cap values ---
  if (scale) mat <- t(scale(t(mat)))
  mat[mat >  cap] <-  cap
  mat[mat < -cap] <- -cap

  # --- Order samples by group ---
  group <- factor(group, levels = unique(group))
  ord_samples <- order(group)
  mat <- mat[, ord_samples, drop = FALSE]
  group_ord <- group[ord_samples]

  # --- Annotations ---
  annotation_row <- data.frame(Group = group_ord)
  rownames(annotation_row) <- colnames(mat)
  quad_vec <- rep(names(quad_genes), lengths(quad_genes))
  annotation_col <- data.frame(Quadrant = factor(quad_vec, levels = names(quad_genes)))
  rownames(annotation_col) <- rownames(mat)

  # --- Colors ---
  if (is.null(palette_group)) {
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      palette_group <- setNames(
        RColorBrewer::brewer.pal(min(8, nlevels(group)), "Set2"),
        levels(group)
      )
    } else {
      palette_group <- setNames(grDevices::hcl.colors(nlevels(group)), levels(group))
    }
  }
  if (is.null(palette_quad)) {
    palette_quad <- c(Q1 = "#2E8B57", Q2 = "#C0392B", Q3 = "#2C3E50", Q4 = "#8E44AD")
  }

  # --- Row/column gaps ---
  row_gaps <- cumsum(table(group_ord))
  col_gaps <- cumsum(lengths(quad_genes))[1:(length(quad_genes)-1)]

  # --- Plot ---
  p <- pheatmap::pheatmap(
    t(mat),
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    annotation_row = annotation_row,
    annotation_col = annotation_col,
    gaps_row = row_gaps,
    gaps_col = col_gaps,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    annotation_colors = list(
      Group    = palette_group,
      Quadrant = palette_quad
    )
  )

  invisible(p)
}
