#' UMAP on samples using selected features
#'
#' @param X numeric matrix [features x samples].
#' @param feats character; feature names to use (non-matching are ignored with a message).
#' @param n_neighbors integer; UMAP neighbor size (default 15).
#' @param min_dist numeric; UMAP min_dist (default 0.1).
#' @param metric character; distance metric (default "euclidean").
#' @param seed integer; RNG seed for reproducibility (default 1).
#' @return matrix [samples x 2] with rownames = sample IDs.
#' @export
run_umap_samples <- function(
    X, feats,
    n_neighbors = 15, min_dist = 0.1,
    metric = "euclidean",
    seed = 1
) {
  stopifnot(is.matrix(X), length(feats) >= 1)
  feats <- intersect(feats, rownames(X))
  if (length(feats) == 0L) stop("None of the requested 'feats' are present in 'X'.")
  S <- scale(t(X[feats, , drop = FALSE]))  # samples x features
  set.seed(seed)
  emb <- uwot::umap(S,
                    n_neighbors = n_neighbors,
                    min_dist = min_dist,
                    metric = metric,
                    ret_model = FALSE)
  rownames(emb) <- rownames(S)
  emb
}


#' Quick ggplot for UMAP embedding
#'
#' @param emb matrix from [run_umap_samples()] (samples x 2).
#' @param groups vector of sample groups (length = nrow(emb)).
#' @param title character; plot title.
#' @param point_size numeric; point size (default 3).
#' @param shape integer; point shape (default 21 = filled circle).
#' @param color,stroke,alpha aesthetics for outlines and transparency.
#' @param palette optional named vector of fill colors; if NULL, generated automatically.
#' @export
plot_umap <- function(
    emb, groups, title = "",
    point_size = 3, shape = 21,
    color = "black", stroke = 1.1, alpha = 0.9,
    palette = NULL
) {
  stopifnot(is.matrix(emb), ncol(emb) >= 2, nrow(emb) == length(groups))
  groups <- factor(groups)
  df <- data.frame(
    UMAP1 = emb[, 1],
    UMAP2 = emb[, 2],
    Group = groups,
    row.names = rownames(emb)
  )

  # Auto palette if not supplied
  if (is.null(palette)) {
    levs <- levels(groups)
    nlev <- length(levs)
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      pal <- RColorBrewer::brewer.pal(min(max(3, nlev), 8), "Set2")
      if (nlev > length(pal))
        pal <- grDevices::rainbow(nlev)
    } else {
      pal <- grDevices::hcl.colors(nlev, "Set2", rev = FALSE)
    }
    palette <- setNames(pal[seq_len(nlev)], levs)
  }

  ggplot2::ggplot(df, ggplot2::aes(UMAP1, UMAP2, fill = Group)) +
    ggplot2::geom_point(size = point_size, alpha = alpha,
                        shape = shape, color = color, stroke = stroke) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(title = title) +
    ggplot2::theme_classic()
}


#' UMAP triptych using delta-ranked features
#'
#' Generates UMAPs for (1) attenuated (+delta), (2) escalated (-delta),
#' and (3) all features. Returns the three embeddings and a combined plot.
#'
#' @param X numeric matrix [features x samples].
#' @param delta named numeric vector of delta scores (names = feature IDs).
#' @param groups factor or character of sample groups (length = ncol(X)).
#' @param top_n integer; number of features to use per side (default 100).
#' @param n_neighbors,min_dist,metric,seed passed to [run_umap_samples()].
#' @param palette optional named color vector for groups; if NULL, generated automatically.
#' @return list with `emb_attenuated`, `emb_escalated`, `emb_all`, and `plot`.
#' @export
umap_triptych_by_delta <- function(
    X, delta, groups,
    top_n = 100,
    n_neighbors = 15, min_dist = 0.1,
    metric = "euclidean",
    seed = 1,
    palette = NULL
) {
  stopifnot(is.matrix(X), length(groups) == ncol(X))
  stopifnot(is.numeric(delta), !is.null(names(delta)))

  feats_common <- intersect(rownames(X), names(delta))
  if (length(feats_common) == 0L)
    stop("No overlap between rownames(X) and names(delta).")

  delta_use <- delta[feats_common]
  ord <- order(delta_use, decreasing = TRUE)
  top_n <- max(1L, min(top_n, length(delta_use)))

  attenuated_feat <- names(delta_use)[ord][seq_len(top_n)]
  escalated_feat  <- names(delta_use)[rev(ord)][seq_len(top_n)]

  emb_attenuated <- run_umap_samples(
    X, feats = attenuated_feat,
    n_neighbors = n_neighbors, min_dist = min_dist,
    metric = metric, seed = seed
  )
  emb_escalated <- run_umap_samples(
    X, feats = escalated_feat,
    n_neighbors = n_neighbors, min_dist = min_dist,
    metric = metric, seed = seed
  )
  emb_all <- run_umap_samples(
    X, feats = feats_common,
    n_neighbors = n_neighbors, min_dist = min_dist,
    metric = metric, seed = seed
  )

  # automatically infer palette if NULL
  if (is.null(palette)) {
    levs <- unique(as.character(groups))
    nlev <- length(levs)
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      pal <- RColorBrewer::brewer.pal(min(max(3, nlev), 8), "Set2")
      if (nlev > length(pal))
        pal <- grDevices::rainbow(nlev)
    } else {
      pal <- grDevices::hcl.colors(nlev, "Set2", rev = FALSE)
    }
    palette <- setNames(pal[seq_len(nlev)], levs)
  }

  p1 <- plot_umap(emb_attenuated, groups, "Samples using attenuated features",
                  palette = palette)
  p2 <- plot_umap(emb_escalated, groups, "Samples using escalated features",
                  palette = palette)
  p3 <- plot_umap(emb_all, groups, "Samples using all features",
                  palette = palette)

  plt <- ggpubr::ggarrange(p1, p2, p3, ncol = 3, nrow = 1,
                           common.legend = TRUE, legend = "right")

  list(
    emb_attenuated = emb_attenuated,
    emb_escalated  = emb_escalated,
    emb_all        = emb_all,
    plot           = plt
  )
}
