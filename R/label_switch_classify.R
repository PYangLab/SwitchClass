# INTERNAL helper (not exported)
.classifyImpt <- function(X, y, t = 25, ntree = 1000, ncores = parallel::detectCores() - 1L) {
  stopifnot(nrow(X) > 1, ncol(X) > 1, length(y) == ncol(X))
  X_t <- t(X)        # samples × features
  y_f <- as.factor(y)
  seeds <- as.integer(sample.int(1e8, t, replace = FALSE))

  run_one <- function(seed_i) {
    set.seed(seed_i)
    rf <- randomForest::randomForest(
      x = X_t,
      y = y_f,
      ntree = ntree,
      importance = TRUE
    )
    rf$importance[, "MeanDecreaseGini"]
  }

  impt_list <- parallel::mclapply(seeds, run_one, mc.cores = max(1L, ncores))
  impt_mat <- do.call(cbind, impt_list)
  rowMeans(impt_mat, na.rm = TRUE)
}

#' Label-switch classification (with optional multicore)
#'
#' Trains two random-forest classifiers with inverted labels and returns
#' \eqn{\delta = I_{rev} - I_{per}} using MeanDecreaseGini importances.
#'
#' @param X numeric matrix [features x samples] with rownames.
#' @param y_reverse character/factor labels (length = ncol(X)) for the reversal scheme.
#' @param y_persist character/factor labels (length = ncol(X)) for the persistence scheme.
#' @param ntimes integer; number of RF repeats to average (default 25).
#' @param seed integer; base seed (default 1).
#' @param ntree integer; trees per forest (default 1000).
#' @param mtry integer or NULL; variables tried at each split (default NULL).
#' @param parallel logical; use multicore (default TRUE).
#' @param ncores integer; cores if `parallel=TRUE` (default detectCores()-1).
#' @return list with named numeric vectors: `delta`, `I_rev`, `I_per`.
#' @export
label_switch_classify <- function(X, y_reverse, y_persist) {
  library(randomForest)
  I_rev <- .classifyImpt(X, y_reverse)
  I_per <- .classifyImpt(X, y_persist)
  list(delta = I_rev - I_per, I_rev = I_rev, I_per = I_per)
}
