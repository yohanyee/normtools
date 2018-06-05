# Construct an empty matrix like another matrix
construct_like_matrix <- function(like_mtx) {
  strucs_source <- as.character(rownames(like_mtx))
  strucs_target <- as.character(colnames(like_mtx))
  mtx <- matrix(NA, nrow=length(strucs_source), ncol=length(strucs_target), dimnames = list(strucs_source, strucs_target))
  return(mtx)
}



# Compute residuals from a given model
residuals_from_model <- function(x, formula, df, train.indices=NULL) {

  # Check that dimensions match
  if (dim(df)[1] != length(x)) {
    stop("df must have same number of observations as length of x")
  }

  # Check train_indices valid and set test_indices
  all.indices <- 1:length(x)
  if (is.null(train.indices)) {
    train.indices <- all.indices
  }
  if (!all(train.indices %in% all.indices)) {
    stop("train.indices are not valid indices")
  }
  test.indices <- setdiff(all.indices, train.indices)

  # Initialize
  if ("RESPONSE_VECTOR_X" %in% colnames(df)) {
    stop("df contains a column with name RESPONSE_VECTOR_X, that conflicts with this code")
  }
  df[["RESPONSE_VECTOR_X"]] <- x
  x_train <- x[train.indices]
  df_train <- df[train.indices,]
  df_test <- df[test.indices,]

  # Fit model
  model_formula <- reformulate(attr(terms(formula), "term.labels"), response="RESPONSE_VECTOR_X")
  model_train <- lm(formula = model_formula, data = df_train )

  # Extract residuals
  residuals_train <- df_train$RESPONSE_VECTOR_X - predict(model_train, newdata=df_train)
  residuals_test <- df_test$RESPONSE_VECTOR_X - predict(model_train, newdata=df_test)

  # Output
  x_out <- numeric(length(x))
  x_out[train.indices] <- residuals_train
  x_out[test.indices] <- residuals_test

  return(x_out)
}
