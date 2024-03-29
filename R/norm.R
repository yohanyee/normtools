#' Successive normalization of rectangular array
#'
#' A method of successively normalizing both rows and columns of a matrix, a la Brad Efron and further described by Olshen et al.
#' See for more details: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2868388/ .
#' Note that the proof on convergence contained within aforementioned is flawed according to the 2013 paper published by the same authors.
#'
#' @param X A rectangular array.
#' @param type Normalization method. This argument is passed as \code{type} to the \code{norm()}.
#' @param tol Tolerance for convergence.
#' @param verbose Be verbose.
#' @param na.set \code{norm()} requires that no elements in the rectangular array be \code{NA}. Replace \code{NA}s with this value.
#' @return Normalized array.
#'
#' @export
norm_successive <- function(X, type="F", tol=1e-8, verbose=TRUE, na.set=0) {
  X <- as.matrix(X)
  X[is.na(X)] <- na.set
  norm_types <- list(f="Frobenius norm", o="One norm", `1`="One norm", i="Infinity norm", m="Maximum modulus", `2`="Spectral (2)-norm")
  norm_difference <- 1
  norm_before <- norm(X, type = type)
  i <- 1
  while(norm_difference >= tol) {
    # Iterate
    rscale <- t(scale(t(X), center = TRUE, scale=TRUE)) # Row scale
    cscale <- scale(rscale, center = TRUE, scale=TRUE) # Column scale

    # Compute norm differences
    norm_after <- norm(cscale, type=type)
    norm_difference <- abs(norm_after - norm_before)
    norm_before <- norm_after
    X <- cscale

    # Print progress
    if (verbose) {
      print(paste("Iteration:", i, "|", norm_types[[tolower(type)]], "difference:", norm_difference))
      i <- i + 1
    }
  }
  return(X)
}



#' Row sum normalization of rectangular array
#'
#' Normalize a rectangular array by dividing each element by the sum of its row.
#' This ensures that each row sum is 1.
#'
#' @param X A rectangular array.
#' @return Normalized array with rows that sum to 1.
#'
#' @export
norm_row_sum <- function(X) {
  return(X/rowSums(X))
}



#' Column sum normalization of rectangular array
#'
#' Normalize a rectangular array by dividing each element by the sum of its column.
#' This ensures that each column sum is 1.
#'
#' @param X A rectangular array.
#' @return Normalized array with columns that sum to 1.
#'
#' @export
norm_col_sum <- function(X) {
  return(X/colSums(X))
}



#' Total sum normalization of rectangular array
#'
#' Normalize a rectangular array by dividing each element by the sum of all elements.
#'
#' @param X A rectangular array.
#' @return Normalized array with all elements summing to 1.
#'
#' @export
norm_full_sum <- function(X) {
  return(X/sum(X))
}



#' Row sum regression of rectangular array
#'
#' Regress out the row sums from elements in a rectangular array.
#' Specifically, each column is replaced by residuals of the linear model that predicts column values by sums of rows.
#' This is similar to dividing by row sums (i.e. \code{norm_rowsum()}) but also takes the differences in slopes (w.r.t. to the row sums) of each column into account.
#'
#' @param X A rectangular array.
#' @return Normalized array of regression residuals.
#'
#' @export
norm_row_regression <- function(X) {
  X_reg <- construct_like_matrix(X)
  X_rsum <- rowSums(X)
  for (i in 1:dim(X)[2]) {
    X_reg[,i] <- residuals(lm(X[,i] ~ X_rsum))
  }
  return(X_reg)
}



#' Column sum regression of rectangular array
#'
#' Regress out the column sums from elements in a rectangular array.
#' Specifically, each row is replaced by residuals of the linear model that predicts row values by sums of columns.
#' This is similar to dividing by column sums (i.e. \code{norm_colsum()}) but also takes the differences in slopes (w.r.t. to the column sums) of each row into account.
#'
#' @param X A rectangular array.
#' @return Normalized array of regression residuals.
#'
#' @export
norm_col_regression <- function(X) {
  X_reg <- construct_like_matrix(X)
  X_csum <- colSums(X)
  for (i in 1:dim(X)[1]) {
    X_reg[i,] <- residuals(lm(X[i,] ~ X_csum))
  }
  return(X_reg)
}


#' Normalization of rectangular array by muiltivariate model fitting
#'
#' Fit a model with all other columns as predictors and extract residuals
#'
#' @param X A rectangular array, or a vector. If a rectangular array is provided, then the model will be applied to each column separately.
#' @param formula a formula describing the model to be fitted, with response variable omitted. If response variable is present, it will be dropped and the appropriate data from \code{X} will be used.
#' @param df data frame with same number of rows as \code{X}, containing the predictors.
#' @param train.indices an optional set of indices corresponding to the rows of \code{X} indicating the data on which the model should be fit. The same model will be then be applied to the remaining indices. If \code{NULL}, then all rows are used in the model.
#' @param progressbar display a progress bar? (default=TRUE)
#' @param cl a cluster object to pass to pbapply if progressbar=TRUE
#' @return Normalized array or vector of regression residuals.
#'
#' @export
norm_by_model <- function(X, formula, df, train.indices=NULL, progressbar=TRUE, cl=NULL) {
  if (is.vector(X)) {
    return(residuals_from_model(X, formula = formula, df=df, train.indices = train.indices))
  } else {
    if (progressbar) {
      if (is.null(cl)) {
        out <- pbapply(X = X, MARGIN = 2, FUN=residuals_from_model, formula=formula, df=df, train.indices=train.indices)
      } else {
        out <- pbapply(X = X, MARGIN = 2, FUN=residuals_from_model, cl=cl, formula=formula, df=df, train.indices=train.indices)
      }
    } else {
      out <- apply(X = X, MARGIN = 2, FUN=residuals_from_model, formula=formula, df=df, train.indices=train.indices)
    }
    return(out)
  }
}
