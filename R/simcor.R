#' Simulates a block correlation matrix
#' 
#' \code{simcor} generates a correlation matrix in the manner of Hardin et al (2013).
#' 
#' @param k Integer. Number of block diagonals.
#' @param size An integer vector of length k containing the number of rows (also columns) of a (square) block
#' diagonal element of the final matrix, which will then be a square matrix of row length \code{sum(size)}.
#' @param rho A numeric vector of length k, containing the mean (off-diagonal) correlations in each block.
#' @param delta Numeric. The mean correlation in the remaining elements of the final matrix, outside the
#' block diagonals.
#' @param epsilon Numeric. The degree of variability in the matrix.
#' @param eidim The shape of the distribution of variability; if \code{eidim==1}, the variability around the
#' mean is exactly epsilon but in random directions. Increasing values of eidim tend to centre the variation.
#' 
#' @references
#' For full details, see
#' \cite{Hardin, J., Garcia, S. R., and Golan, D. (2013). A method for generating realistic correlation matrices. Annals of Applied Statistics, 7(3):1733-1762.}
#' This function is not intended for use by the user.
#' 
simcor <- function(k = 6, size = c(10,5,8,2,15,50), rho = c(0.7,0.7,0.5,0.9,0.85,0.4),
                   delta = 0.39, epsilon = 0.99 - max(rho), eidim = 2){
  ndim <- sum(size)
  bigcor <- matrix(rep(delta, ndim * ndim), ncol = ndim)
  for (i in 1:k) {
    cor <- matrix(rep(rho[i], size[i] * size[i]), ncol = size[i])
    if (i == 1) {bigcor[1:size[1], 1:size[1]] <- cor}
    if (i != 1) {bigcor[(sum(size[1:(i - 1)]) + 1):sum(size[1:i]),
                        (sum(size[1:(i - 1)]) + 1):sum(size[1:i])] <- cor}
  }
  diag(bigcor) <- 1 - epsilon
  eivect <- c()
  for (i in 1:ndim) {
    ei <- runif(eidim, -1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei/sqrt(sum(ei^2)))
  }
  bigE <- t(eivect) %*% eivect
  cor.nz <- bigcor + bigE
  cor.nz
}