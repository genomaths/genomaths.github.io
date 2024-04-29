#' Modulo Operation
#' @rdname mod
#' @aliases mod
#' @aliases '%%'
#' @aliases modulo
#' @param n A numeric vector (preferably of integers), a matrix where each 
#' element can be reduced to integers.
#' @param m An integer vector (positive, zero, or negative).
#' @description 
#' Integer remainder of the division of the integer n by m:
#' n mod m. 
#' @param ... Not in use.
#' @return An element of x, an \code{\link{Automorphism-class}} object.
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @examples 
#' ## Example 1
#' ## Build a matrix 'n' and set a vector of integers 'm'
#' n <- diag(x=1, nrow = 4, ncol = 4) * c(43,125,2,112)
#' m <- c(64,4,4,64)
#' 
#' ## Operation n mod m 
#' mod(n = n, m = m)
#' 
#' ## Or simply:
#' n %% m
#' 
#' ## Example 2
#' m <- matrix(c(8,2,3, 11,12,13), nrow = 2)
#' m
#' 
#' m %% 4
#' 
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







