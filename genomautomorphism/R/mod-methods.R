#' Modulo Operation
#' @rdname mod
#' @aliases mod
#' @aliases '%%'
#' @aliases modulo
#' @param n A matrix where each element can be reduced to integers or the
#' same as in \code{\link[mod]{numbers}}.
#' @param m As in \code{\link[mod]{numbers}}.
#' @param ... Not in use yet.
#' @description Integer remainder of the division of the integer n by m:
#' n mod m. This function extend the application of function 
#' \code{\link[mod]{numbers}} to matrices where the operation on each row is 
#' with is accomplish with a different values of \eqn{m}, i.e, where \eqn{m}
#' is a vector. 
#' @return An element of x, an \code{\link{Automorphism-class}} object.
#' @importFrom numbers mod
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @examples 
#' ## Build a matrix 'n' and set a vector of integers 'm'
#' n <- diag(x=1, nrow = 4, ncol = 4) * c(43,125,2,112)
#' m <- c(64,4,4,64)
#' 
#' ## Operation n mod m 
#' mod(n = n, m = m)
#' 
#' ## Or simply:
#' n %% m
setGeneric(
    "mod",
    function(n, m, ...) {
        n %% m
    }
)

#' @rdname mod
#' @aliases mod
#' @aliases '%%'
#' @aliases modulo
setMethod(
    "mod", signature(n = "matrix", m = "numeric"),
    function(n, m) {
        m <- as.integer(m)
        if (length(m) < 2) 
            res <- n %% m
        else
            res <- mapply(function(x, y) x %% y, 
                split(n, row(n)), m, USE.NAMES = FALSE)
        return(res)
    }
)







