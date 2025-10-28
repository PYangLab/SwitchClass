#' Assign SwitchClass quadrants
#'
#' Q1:  delta >  T & FC >  0
#' Q2:  delta < -T & FC >  0
#' Q3:  delta < -T & FC <  0
#' Q4:  delta >  T & FC <  0
#' Q0:  otherwise
#'
#' @param delta named numeric vector (features).
#' @param fc    named numeric vector (same feature names as `delta`).
#' @param delta_thresh numeric threshold for |delta|.
#' @return factor vector of quadrant labels with names = features.
#' @export
assign_quadrants <- function(delta, fc, delta_thresh = 0.15) {
  common <- intersect(names(delta), names(fc))
  if (length(common) == 0L) stop("No overlapping names between 'delta' and 'fc'.")
  d <- delta[common]; f <- fc[common]

  q <- ifelse(d >  delta_thresh & f >  0, "Q1",
              ifelse(d < -delta_thresh & f >  0, "Q2",
                     ifelse(d < -delta_thresh & f <  0, "Q3",
                            ifelse(d >  delta_thresh & f <  0, "Q4", "Q0"))))
  structure(factor(q, levels = c("Q1","Q2","Q3","Q4","Q0")), names = common)
}


#' Scatter plot of delta vs fold-change with quadrant assignment
#'
#' Builds a tidy data.frame with quadrants and returns the ggplot scatter.
#'
#' @param delta named numeric vector of δ (features).
#' @param fc    named numeric vector of log2 fold-change (same names as `delta`).
#' @param delta_thresh numeric; |δ| threshold to define Q1–Q4 (default 0.15).
#' @param label_top integer; number of labels per quadrant (by |δ|). Set 0 to skip.
#' @param label_size numeric; text size for labels.
#' @param palette named colors for quadrants. Names should include Q1..Q4 and optional Q0.
#' @param point_size numeric; scatter point size.
#' @param alpha numeric in [0,1]; point transparency.
#' @param xlab,ylab, title axis labels and plot title.
#' @return list with elements: `data` (data.frame) and `plot` (ggplot object).
#' @export
plot_delta_vs_fc <- function(
    delta, fc,
    delta_thresh = 0.15,
    label_top = 5,
    label_size = 3,
    palette = c(Q1 = "#2E8B57", Q2 = "#C0392B", Q3 = "#C0392B", Q4 = "#2E8B57", Q0 = "grey70"),
    point_size = 2.2,
    alpha = 0.85,
    xlab = expression(delta),
    ylab = expression(log[2]*" FC"),
    title = "feature importance vs log2 fold-change"
) {
  stopifnot(is.numeric(delta), is.numeric(fc))
  quads <- assign_quadrants(delta, fc, delta_thresh = delta_thresh)
  feats <- names(quads)
  df <- data.frame(
    feature = feats,
    delta   = as.numeric(delta[feats]),
    fc      = as.numeric(fc[feats]),
    quadrant = factor(as.character(quads), levels = levels(quads)),
    stringsAsFactors = FALSE
  )

  # choose labels: top |δ| per quadrant (skip Q0)
  label_df <- NULL
  if (label_top > 0) {
    label_df <- do.call(rbind, lapply(c("Q1","Q2","Q3","Q4"), function(q) {
      sub <- df[df$quadrant == q & is.finite(df$delta) & is.finite(df$fc), , drop = FALSE]
      if (nrow(sub) == 0) return(NULL)
      sub[order(abs(sub$delta), decreasing = TRUE), ][seq_len(min(label_top, nrow(sub))), ]
    }))
  }

  # ensure palette covers observed quadrants
  present_quads <- levels(df$quadrant)[levels(df$quadrant) %in% unique(as.character(df$quadrant))]
  pal_use <- palette
  missing_cols <- setdiff(present_quads, names(pal_use))
  if (length(missing_cols)) {
    pal_use <- c(pal_use, setNames(rep("grey60", length(missing_cols)), missing_cols))
  }

  # plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = delta, y = fc, color = quadrant)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", size = 0.3, color = "grey60") +
    ggplot2::geom_vline(xintercept =  c(-delta_thresh, delta_thresh), linetype = "dotted", size = 0.3, color = "grey60") +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::scale_color_manual(values = pal_use, drop = FALSE) +
    ggplot2::labs(x = xlab, y = ylab, title = title, color = "Quadrant") +
    ggplot2::theme_classic()

  if (!is.null(label_df) && nrow(label_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = feature),
      size = label_size, max.overlaps = Inf, show.legend = FALSE
    )
  }

  list(data = df, plot = p)
}
