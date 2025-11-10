RUV4_upgrade <- function(Y, X, ctl, k, Z = NULL, eta = NULL, include.intercept = TRUE, 
                 fullW0 = NULL, inputcheck = TRUE) 
{
  m = nrow(Y)
  n = ncol(Y)
  ctl = ctl2logi(ctl, n)
  Y = ruv::RUV1(Y, eta, ctl, include.intercept = include.intercept)
  # If covariates Z are provided, residualize Y with respect to both X and Z
  if (!is.null(Z)) {
    # Ensure Z is a matrix with m rows
    if (is.vector(Z)) {
      Z <- matrix(Z, ncol = 1)
    }
    if (!is.matrix(Z)) {
      Z <- as.matrix(Z)
    }
    if (nrow(Z) != m) {
      stop("RUV4_upgrade: Number of rows in Z must match number of samples (rows of Y).")
    }
    XZ <- cbind(X, Z)
    Y0 = orthogonal_projection(Y, XZ)
    p_eff <- ncol(XZ)
  } else {
    Y0 = orthogonal_projection(Y, X)
    p_eff <- ncol(X)
  }
  if (is.null(fullW0)) {
    full_U = svd(Y0 %*% t(Y0))$u[, 1:(m - p_eff), drop = FALSE]
  }
  if (k > 0) {
    U = full_U[, 1:k, drop = FALSE]
    alpha = t(U) %*% Y0
    Y0c = Y0[, ctl, drop = FALSE]
    W = Y0c %*% t(Y0c) %*% U %*% solve(t(U) %*% Y0c %*% t(Y0c) %*% U)
  }
  # Remove unwanted variation captured by W
  newY = Y - W %*% alpha
  # If covariates provided, also regress out Z from the corrected matrix
  if (!is.null(Z)) {
    # Ordinary least squares projection to remove Z effects
    newY = orthogonal_projection(newY, Z)
  }
  return(list(
    newY = newY,
    full_U = full_U,
    W = W,
    WA = W %*% alpha,
    alpha = alpha))
}

ctl2logi <- function(ctl, n)
{
  ctl2 = rep(FALSE, n)
  ctl2[ctl] = TRUE
  return(ctl2)
}

orthogonal_projection <- function (A, B) {
  return(A - B %*% solve(t(B) %*% B) %*% t(B) %*% A)
}