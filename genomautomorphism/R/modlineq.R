## Copyright (C) 2021 Robersy Sanchez <https://genomaths.com/>
## Author: Robersy Sanchez This file is part of the R package
## 'GenomAutomorphism'.  'GenomAutomorphism' is a free
## software: you can redistribute it and/or modify it under the
## terms of the GNU General Public License as published by the Free
## Software Foundation, either version 3 of the License, or (at
## your option) any later version.  This program is distributed in
## the hope that it will be useful, but WITHOUT ANY WARRANTY;
## without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for
## more details.  You should have received a copy of the GNU
## General Public License along with this program; if not, see
## <http://www.gnu.org/licenses/>.

#' @rdname modlineq
#' @aliases modlineq
#' @title Modular System of Linear Equation Solver (MLE)
#' @param a An integer or a vector of integers. 
#' @param b An integer or a vector of integers.
#' @param n An integer or a vector of integers.
#' @param no.sol Values to return when the equation is not solvable or yield
#' the value 0. Default is 0.
#' @description If \eqn{a, b}, and  \eqn{c} are integer vectors, this function 
#' try to find, at each coordinate, the solution of the MLE 
#' \eqn{a x = b}  mod \eqn{n}. If the MLE \eqn{a x = b mod n} has not 
#' solutions (see \code{\link[numbers]{modlin}}), the value reported for the 
#' coordinate will be 0 and the corresponding translation.
#' @details For \eqn{a, b}, and \eqn{c} integer scalars, it is just a 
#' wrapper function to call \code{\link[numbers]{modlin}}. 
#' @importFrom numbers modlin
#' @export
#' @return If the solution is exact, then a numerical vector will be returned,
#' otherwise, if there is not exact solution for some coordinate, the a list
#' carrying the element on the diagonal matrix and a translation vector will
#' be returned.
#' @examples
#' ## Set the vector x, y, and m.
#' x <- c(9,32,24,56,60,27,28,5)
#' y <-  c(8,1,0,56,60,0,28,2)
#' modulo <- c(64,125,64,64,64,64,64,64)
#' 
#' ## Try to solve the modular equation a x = b mod n
#' m <- modlineq(a = x, b = y, n = modulo)
#' m
#' 
#' ## Or in matrix form 
#' diag(m)
#' 
#' ## The reverse mapping is an affine transformation
#' mt <- modlineq(a = y, b = x, n = modulo, no.sol = 1L)
#' mt
#' 
#' ## That is, vector 'x' is revovered with the transformaiton
#' (y %*% diag(mt$diag) + mt$translation) %% modulo
#' 
#' # Or
#' cat("\n---- \n")
#' 
#' (y %*% diag(mt$diag) + mt$translation) %% modulo == x
setGeneric(
    "modlineq",
    function(a, b, n, no.sol = 0L) {

        if (is.numeric(no.sol))
            no.sol <- as.integer(no.sol)
        
        if (length(a) < 2) 
            res <- modl(a, b, n, no.sol = no.sol)
        else {
            if (length(n) == 1) {
                res <- mapply(function(x,y) modl(x, y, n, no.sol = no.sol),
                            a, b, USE.NAMES = FALSE)
            }
            else
                res <- mapply(function(x, y, z) {
                    modl(x, y, z, no.sol = no.sol)
                }, a, b, n, USE.NAMES = FALSE)
        }
        
        ## ========= Checking ===========
        trls <- FALSE
        idx <- (a %*% diag(res)) %% n != b
        if (sum(idx, na.rm = TRUE) > 0) {
            trls <- TRUE
            ## Compute translations
            if (is.integer(no.sol))
                trl <- c(b - (a %*% diag(res)) %% n) %% n
        }
        
        if (trls) 
            res <- list(diag = res, translation = trl)
            
        return(res)
    }
)

## ========= Auxiliary function ===================

modl <- function(x, y, m, no.sol = 0L) {
    x <- as.integer(x)
    y <- as.integer(y)
    m <- as.integer(m)
    
    if (length(m) == 1 && length(x) > 1)
        m <- rep(m, length(x))
    if (x > 0 && y > 0) {
        res <- modlin(x, y, m)
        if (length(res) > 1) 
            res <- min(res)
    }
    else 
        res <- no.sol
    if (is.null(res)) 
        res <- no.sol
    
    return(res)
}



