#' Add noise to any user specified correlation matrix
#' 
#' \code{noisecor} adds noise to any user specified correlation matrix, in the manner of Hardin et al (2013).
#' 
#' @param cormat A correlation matrix. The function does not itself check whether this is a valid
#' correlation matrix.
#' @param epsilon Numeric. The maximim perturbation of any correlation.
#' @param eidim Integer. The shape of the perturbation; if \code{eidim==1}, the perturbation is exactly
#' epsilon but in random directions. Increasing values of eidim tend to centre the perturbations.
#' 
#' @references 
#' For full details, see
#' \cite{Hardin, J., Garcia, S. R., and Golan, D. (2013). A method for generating realistic correlation matrices. Annals of Applied Statistics, 7(3):1733-1762.}
#' This function is not intended for use by the user.
#'
noisecor <- function(cormat, epsilon = .01, eidim=2){
  ndim <- dim(cormat)[1]
  diag(cormat) <- 1 - epsilon
  eivect <- c( )
  for (i in 1:ndim) {
    ei <- runif(eidim, -1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei^2)) )
  }
  bigE <- t(eivect) %*% eivect
  cor.nz <- cormat + bigE
  cor.nz
}