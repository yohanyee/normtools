% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm.R
\name{norm_row_regression}
\alias{norm_row_regression}
\title{Row sum regression of rectangular array}
\usage{
norm_row_regression(X)
}
\arguments{
\item{X}{A rectangular array.}
}
\value{
Normalized array of regression residuals.
}
\description{
Regress out the row sums from elements in a rectangular array.
Specifically, each column is replaced by residuals of the linear model that predicts column values by sums of rows.
This is similar to dividing by row sums (i.e. \code{norm_rowsum()}) but also takes the differences in slopes (w.r.t. to the row sums) of each column into account.
}
