#' Reactome ORA per quadrant (using pathwayOverrepresent)
#'
#' Wraps a convenience pipeline:
#' - builds Reactome TERM2GENE via msigdbr
#' - runs pathwayOverrepresent on each quadrant (Q1..Q4)
#' - returns top-k barplots (−log10 p) and tidy tables
#'
#' @param df data.frame with columns: gene, delta, quadrant (Q1..Q4)
#' @param universe character vector of background gene symbols (e.g., rownames(discovery.cohort))
#' @param species msigdbr species (default "Homo sapiens")
#' @param top_terms how many top terms per quadrant to plot (default 5)
#' @param palette named colors for bars of Q1..Q4
#' @return list with:
#'   \describe{
#'     \item{results}{named list of data.frames p1..p4 (or NULL if no hits)}
#'     \item{barplot_df}{tidy data.frame used for plotting}
#'     \item{plot}{ggpubr patch of 4 panels}
#'   }
#' @export
enrich_reactome_quadrants <- function(
    df, universe,
    species = "Homo sapiens",
    top_terms = 5,
    palette = c(Q1 = "#33DBB2", Q2 = "#DFB647", Q3 = "#2ECEF8", Q4 = "#76CF49")
) {

  df$quadrant <- ifelse(df$delta >  0 & df$fc >  0, "Q1",
                        ifelse(df$delta < 0 & df$fc >  0, "Q2",
                               ifelse(df$delta < 0 & df$fc <  0, "Q3",
                                      ifelse(df$delta >  0 & df$fc <  0, "Q4", "Q0"))))

  stopifnot(all(c("feature","quadrant") %in% names(df)))
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Package 'msigdbr' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Package 'dplyr' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")
  if (!requireNamespace("ggpubr", quietly = TRUE))
    stop("Package 'ggpubr' is required.")

  # Build Reactome TERM2GENE via msigdbr
  msig_react <- msigdbr::msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME") |>
    dplyr::select(gs_name, gene_symbol)
  rect <- split(msig_react$gene_symbol, msig_react$gs_name)

  # Helper to run ORA via pathwayOverrepresent safely
  run_overrep <- function(genes, annotation, universe) {
    genes <- intersect(genes, unique(unlist(annotation)))
    if (length(genes) == 0) return(NULL)
    pathwayOverrepresent(genes, annotation = annotation, universe = universe, alter = "greater")
  }

  # Per quadrant gene sets (names are symbols already)
  quads <- c("Q1","Q2","Q3","Q4")
  g_by_q <- lapply(quads, function(q) df$feature[df$quadrant == q])
  names(g_by_q) <- quads

  # Run ORA for each quadrant
  p_list <- lapply(g_by_q, run_overrep, annotation = rect, universe = universe)
  names(p_list) <- quads   # p_list$Q1 etc.

  # Build the plotting df (handle missing/short results)
  top_pick <- function(p_tab, k) {
    k <- min(k, nrow(p_tab))
    out <- as.data.frame(p_tab)
    out$Term <- rownames(out)
    out$pvalue <- as.numeric(out$pvalue)
    out <- out[order(out$pvalue, decreasing = FALSE), , drop = FALSE]
    out <- out[seq_len(k), , drop = FALSE]
    return(out)
  }

  p1 <- top_pick(p_list$Q1, top_terms)
  p2 <- top_pick(p_list$Q2, top_terms)
  p3 <- top_pick(p_list$Q3, top_terms)
  p4 <- top_pick(p_list$Q4, top_terms)

  # Tidy barplot_df
  as_bar_df <- function(ptab, label) {
    if (is.null(ptab)) return(data.frame(Term = character(0), Value = numeric(0), Group = character(0)))
    data.frame(
      Term  = ptab$Term,
      Value = -log10(as.numeric(ptab[,"pvalue"])),
      Group = label,
      stringsAsFactors = FALSE
    )
  }

  barplot_df <- dplyr::bind_rows(
    as_bar_df(p1, "Q1"),
    as_bar_df(p2, "Q2"),
    as_bar_df(p3, "Q3"),
    as_bar_df(p4, "Q4")
  ) %>%
    dplyr::mutate(
      Term = gsub("^REACTOME[_:]*", " ", Term)
    )

  # Make 4 bar plots (empty panels if no terms)
  make_one <- function(df_sub, qlab) {
    ggplot2::ggplot(df_sub, ggplot2::aes(x = stats::reorder(Term, -Value), y = Value)) +
      ggplot2::geom_bar(stat = "identity", fill = palette[[qlab]]) +
      ggplot2::labs(title = qlab, x = NULL, y = expression(-log[10]*" p")) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        plot.background  = ggplot2::element_blank(),
        panel.border     = ggplot2::element_blank(),
        panel.grid       = ggplot2::element_blank(),
        axis.text.x      = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
  }

  pQ1 <- make_one(dplyr::filter(barplot_df, Group == "Q1"), "Q1")
  pQ2 <- make_one(dplyr::filter(barplot_df, Group == "Q2"), "Q2")
  pQ3 <- make_one(dplyr::filter(barplot_df, Group == "Q3"), "Q3")
  pQ4 <- make_one(dplyr::filter(barplot_df, Group == "Q4"), "Q4")

  combined <- ggpubr::ggarrange(pQ1, pQ2, pQ3, pQ4, ncol = 4, nrow = 1)

  list(
    results = p_list,         # raw p1..p4 tables (from pathwayOverrepresent)
    barplot_df = barplot_df,  # tidy, combined top-terms table
    plot = combined           # 4-panel barplot
  )
}
