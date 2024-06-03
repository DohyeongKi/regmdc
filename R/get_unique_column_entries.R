#' Get the unique column entries of a design matrix
#'
#' Given a design matrix, this function returns the list whose ith element is
#' the sorted vector of the unique entries appearing in the ith column of the
#' matrix. For entirely monotonic regression, Hardy—Krause variation denoising,
#' and their generalization ('emhk'), 0 is dropped if exists. For totally
#' concave regression, MARS via LASSO, and their generalization ('tcmars'), 0 is
#' additionally included and the maximal element is dropped.
#'
#' @param X_design A numeric design matrix. Each row corresponds to individual
#'   data.
#' @param method A string indicating the estimation method and determining the
#'   post-processing. One of "emhk" and "tcmars".
#' @param number_of_bins An integer vector of the numbers of bins for the
#'   approximate methods. `NULL` if the approximate methods are not used.
#'   Currently available for "tc", "mars", and "tcmars".
#'
#' @details
#' If `number_of_bins` is not `NULL`, then instead of the exact values of
#' entries, the endpoints of the equally-spaced bins constructed according to
#' `number_of_bins` are exploited.
get_unique_column_entries <- function(X_design, method, number_of_bins = NULL) {
  if (method %in% c('emhk')) {
    lapply((1L:ncol(X_design)), function(col) {
      column_unique <- unique(X_design[, col])
      column_unique <- column_unique[column_unique > 0]
      sort(column_unique)
    })
  } else if (method %in% c('tcmars')) {
    if (is.null(number_of_bins)) {
      lapply((1L:ncol(X_design)), function(col){
        column_unique <- unique(c(0, X_design[, col]))
        sort(column_unique)
        column_unique <- column_unique[-length(column_unique)]
      })
    } else {
      lapply((1L:ncol(X_design)), function(col){
        N <- number_of_bins[col]
        if (is.na(N)) {
          column_unique <- unique(c(0, X_design[, col]))
          sort(column_unique)
          column_unique <- column_unique[-length(column_unique)]
        } else {
          column_unique <- seq(0L, N - 1L) / N
        }
      })
    }
  } else {
    stop('`method` must be one of "emhk" and "tcmars".')
  }
}
