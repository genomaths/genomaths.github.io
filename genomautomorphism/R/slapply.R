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

#' @rdname slapply
#' @aliases slapply
#' @title Apply a function over a list-like object preserving its attributes
#' @param x A list-like or vector-like object.
#' @param FUN,... The same as described in \code{\link[base]{lapply}}.
#' @param keep.attr Logic. If TRUE, then the original attributes from 'x' are
#' preserved in the returned list. Default is FALSE.
#' @param class Name of the class to which the returned list belongs to.
#' Default is NULL.
#' @param simplify,USE.NAMES The same as described in
#' \code{\link[base]{sapply}}.
#' @description This function apply a function over a list-like object
#' preserving its attributes and simplify (if requested) the list as
#' \code{\link[base]{sapply}} function does. \strong{slapply} returns a
#' list of the same length as 'x', each element of which is the result of
#' applying FUN to the corresponding element of 'x'.
#' @return Same as in ?base::\code{\link[base]{lapply}} if keep.attr = FALSE.
#' Otherwise same values preserving original attributes from 'x'.
#' @seealso \code{\link[base]{lapply}} and \code{\link[base]{sapply}}
#' @export
#' @examples
#' ## Create a list
#' x <- list(a = seq(10), beta = exp(seq(-3, 3)), 
#'           logic = c(TRUE, FALSE, FALSE, TRUE))
#' class(x) <- "nice"
#'
#' ## To compute the list mean for each list element using 'base::lapply'
#' class(slapply(x, mean, simplify = FALSE))
#'
#' ## Simply 'base::lapply' preserving attributes
#' slapply(x, mean, keep.attr = TRUE, simplify = FALSE)
#'
#' ## To preserve attributes and simplify
#' slapply(x, mean, keep.attr = TRUE, simplify = TRUE)
#'
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @export
#' @export
slapply <- function(x,
    FUN,
    keep.attr = FALSE,
    class = NULL,
    simplify = TRUE,
    USE.NAMES = TRUE,
    ...) {
    if (keep.attr) {
        s4 <- typeof(x) == "S4"
        nm <- names(x)
        if (is.null(class)) {
            cl <- class(x)
        }
        x <- lapply(x, FUN, ...)

        if (simplify) {
            x <- unlist(x)
        }
        if (!USE.NAMES) {
            x <- unname(x)
        }

        if (s4) {
            names(x) <- nm
            if (!is.null(class)) {
                x <- as(x, class)
            } else {
                y <- try(as(x, cl))
                if (!inherits(x, "try-error")) {
                    x <- y
                }
            }
        } else {
            if (!is.null(class)) {
                x <- structure(x, class = class, names = nm)
            } else {
                y <- try(structure(x, class = cl, names = nm))
                if (!inherits(y, "try-error")) {
                    x <- y
                }
            }
        }
    } else {
        x <- lapply(x, FUN, ...)
        if (simplify) {
            x <- unlist(x)
        }
        if (!USE.NAMES) {
            x <- unname(x)
        }
    }
    return(x)
}
