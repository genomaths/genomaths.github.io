## Copyright (C) 2021-2024 Robersy Sanchez <https://genomaths.com/>
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

## ======================='[' AutomorphismList ======================


#' @rdname extract-methods
#' @aliases '['
#' @aliases extract
#' @aliases extract-methods
#' @title An S4 class to extract elements for objects created with
#' GenomAutomorphism package
#' @param x An object from [AutomorphismList], [ListCodonMatrix], or 
#' [MatrixSeq].
#' @param i,j,... As in \code{\link[base]{Extract}}.
#' @description First and second level subsetting of 'x'. Extraction using
#' names can be done as x$name.
#' @return An object from [AutomorphismList], [ListCodonMatrix], or 
#' [MatrixSeq] class.
#' @keywords internal
#' @exportMethod "["
#' @export
#' @author Robersy Sanchez <https://genomaths.com>
setMethod(
    "[", signature(x = "AutomorphismList"),
    function(x, i, ...) {
        x@DataList <- x@DataList[i]
        return(x)
    }
)

## ======================='[' ListCodonMatrix ======================

#' @rdname extract-methods
#' @aliases '['
#' @aliases extract
#' @aliases extract-methods
#' @keywords internal
#' @exportMethod "["
#' @export
#' @author Robersy Sanchez <https://genomaths.com>
setMethod(
    "[", signature(x = "ListCodonMatrix"),
    function(x, i, ...) {
        x@DataList <- x@DataList[i]
        x@seq_alias <- x@seq_alias[i]
        return(x)
    }
)

## ======================= '[' MatrixSeq ======================

#' @rdname extract-methods
#' @aliases '['
#' @aliases extract
#' @aliases extract-methods
#' @keywords internal
#' @exportMethod "["
#' @export
#' @author Robersy Sanchez <https://genomaths.com>
setMethod(
    "[", signature(x = "MatrixSeq"),
    function(x, i, j, ...) {
        cn <- colnames(x@matrix)
        rn <- rownames(x@matrix)
        if (missing(j)) {
            x@matrix <- as.matrix(x@matrix[i, ])
            if (length(i) == 1)
                colnames(x@matrix) <- rn[i]
        }
        if (missing(i)) {
            x@matrix <- as.matrix(x@matrix[, j])
            colnames(x@matrix) <- cn[j]
        }
        
        if (!missing(i) && !missing(j)) {
            x@matrix <- as.matrix(x@matrix[i, j])
            colnames(x@matrix) <- cn[j]
            rownames(x@matrix) <- rn[i]
        }
        x@seqs <- x@seqs[i]
        x@names <- x@names[i]
        return(x)
    }
)


## ======================= '[[' AutomorphismList ======================


#' @rdname extract-methods
#' @aliases '[['
#' @aliases extract
#' @aliases extract-methods
#' @exportMethod "[["
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @examples
#' ## Load automorphisms found BRCA1 primate genes
#' data("brca1_autm", package = "GenomAutomorphism")
#' 
#' ## Extract AutomorphismList object with only one element
#' brca1_autm[1]
#' 
#' ## Extract Automorphism object with only one element
#' brca1_autm[[3]]
#' 
#' ## Extract Automorphism object using element name.
#' brca1_autm[["human_1.gorilla_1"]]
setMethod(
    "[[", signature(x = "AutomorphismList"),
    function(x, i, ...) {
        x <- x[i]
        x <- getAutomorphisms(x)
        x <- as(x, "GRangesList")
        return(x[[1]])
    }
)


## ===================== '[[' ListCodonMatrix ======================


#' @rdname extract-methods
#' @aliases '[['
#' @aliases extract
#' @aliases extract-methods
#' @exportMethod "[["
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
setMethod(
    "[[", signature(x = "ListCodonMatrix"),
    function(x, i, ...) {
        x <- x[i]
        x <- x@DataList[[1]] 
        return(x)
    }
)

## ======================= $ AutomorphismList ======================

#' @rdname extract-methods
#' @aliases extract
#' @aliases $
#' @aliases extract-methods
#' @param name Element name in the list 'x'.
#' @exportMethod "$"
#' @export
setMethod(
    "$", signature(x = "AutomorphismList"),
    function(x, name) {
        x@DataList <- x@DataList[match(name, names(x))]
        x <- getAutomorphisms(x)
        x <- as(x, "GRangesList")
        return(x[[1]])
    }
)


## ======================= $ ListCodonMatrix ======================

#' @rdname extract-methods
#' @aliases extract
#' @aliases extract-methods
#' @exportMethod "$"
#' @export
setMethod(
    "names", signature(x = "ListCodonMatrix"),
    function(x) {
        x@names
    }
)

#' @rdname extract-methods
#' @aliases $
#' @aliases extract
#' @aliases extract-methods
#' @exportMethod "$"
#' @export
setMethod(
    "$", signature(x = "ListCodonMatrix"),
    function(x, name) {
        i <- match(name, names(x))
        x <- x[i]
        x <- x@DataList[[1]] 
        return(x)
    }
)


## ======================= $ MatrixSeq ======================


#' @rdname extract-methods
#' @aliases $
#' @aliases extract
#' @aliases extract-methods
#' @exportMethod "$"
#' @export
setMethod(
    "$", signature(x = "MatrixSeq"),
    function(x, name) {
        i <- match(name, names(x))
        return(x[i])
    }
)


#' @rdname extract-methods
#' @aliases extract
#' @aliases extract-methods
#' @param value A character vector of up to the same length as x, or NULL.
#' @exportMethod "$"
#' @export
setMethod(
    "names<-", signature(x = "MatrixSeq"),
    function(x, value) {
        if (length(x@names) == length(value)) {
            x@names <- value
            if (is.matrix(x@matrix))
                rownames(x@matrix) <- value
            else
                names(x@matrix) <- value
            
            names(x@seqs) <- value
        }
        return(x)
    }
)















