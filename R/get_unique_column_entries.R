#' Get the unique column entries of a design matrix
#'
#' Given a design matrix, this function returns the list whose ith element is
#' the sorted vector of the unique entries appearing in the ith column of the
#' matrix. For entirely monotonic regression, Hardy—Krause variation denoising,
#' and their generalization ('emhk'), 0 is dropped if exists. For totally convex
#' regression, MARS via LASSO, and their generalization ('tcmars'), 0 is
#' additionally included and the maximal element is dropped.
#'
#' @param X_design A numeric design matrix. Each row corresponds to individual
#'   data.
#' @param method A string indicating the estimation method and determining the
#'   post-processing. One of "emhk" and "tcmars".
get_unique_column_entries <- function(X_design, method) {
  if (method %in% c('emhk')) {
    lapply((1L:ncol(X_design)), function(col) {
      column_unique <- unique(X_design[, col])
      column_unique <- column_unique[column_unique > 0]
      sort(column_unique)
    })
  } else if (method %in% c('tcmars')) {
    lapply((1L:ncol(X_design)), function(col){
      column_unique <- unique(c(0, X_design[, col]))
      sort(column_unique)
      column_unique <- column_unique[-length(column_unique)]
    })
  } else {
    stop('`method` must be one of "emhk" and "tcmars".')
  }
}
