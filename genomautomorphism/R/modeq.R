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

#' @rdname modeq
#' @aliases modeq
#' @title A Wrapper Calling Modular Linear Equation Solver (MLE)
#' @description It is just a wrapper function to call
#' \code{\link[numbers]{modlin}}. This function is intended to be use
#' internally. MLE (\eqn{a * x = b mod n}) not always has solution If the MLE
#' has not solution the function will return the value -1. Also, if
#' \eqn{a * x = b mod n} has solution x = 0, then function \emph{'modeq'}
#' will return -1.
#' @keywords internal
#' @importFrom numbers modlin
#' @export
#' @return A number. If the equation has not solution in their definition, 
#' domain it will return -1.
#' @examples
#' ## The MLE 10 * x = 3 mod 64 has not solution
#' modeq(10, 3, 64)
#'
#' ## The result is the giving calling modlin(10, 4, 64)
#' modeq(10, 4, 64)
modeq <- function(a, b, n) {
    if (all(slapply(c(a, b, n), is.numeric))) {
        a <- as.integer(a)
        b <- as.integer(b)
        n <- as.integer(n)
        
        if (!any(is.na(c(a, b)))) {
            if (a > 0 && b > 0) {
                res <- modlin(a, b, n)[1] ## first solution
            } else {
                if (a != b)
                    res <- -1
            }
        } else {
            res <- -1
        }
    } else {
        res <- -1
    }

    if (is.null(res) || res == 0) {
        res <- -1
    }
    return(res)
}

