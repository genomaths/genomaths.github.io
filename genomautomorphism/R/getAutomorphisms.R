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

#' @rdname getAutomorphisms
#' @aliases getAutomorphisms
#' @title Get Automorphisms
#' @param x An \code{\link{AutomorphismList-class}}.
#' @param ... Not in use.
#' @description For the sake of saving memory, each
#' \code{\link{Automorphism-class}}
#' objects is stored in an \code{\link{AutomorphismList-class}}, which  does
#' not inherits from a \code{\link[GenomicRanges]{GRanges-class}}.
#' @details This function just transform each \code{\link{Automorphism-class}}
#' object into an object from the same class but now inheriting from a
#' \code{\link[GenomicRanges]{GRanges-class}}.
#' @export
#' @examples
#' ## Load a dataset
#' data(autm, package = "GenomAutomorphism")
#' aut <- mcols(autm)
#' aut ## This a DataFrame object
#'
#' ## The natural ranges for the sequence (from 1 to length(aut)) are added
#' getAutomorphisms(aut)
#'
#' ## A list of automorphisms
#' aut <- list(aut, aut)
#' getAutomorphisms(aut)
#'
#' ## Automorphism-class inherits from 'GRanges-class'
#' aut <- as(autm, "GRanges")
#' as(aut, "Automorphism")
setGeneric(
    "getAutomorphisms",
    function(x,
    ...) {
        standardGeneric("getAutomorphisms")
    }
)


#' @rdname getAutomorphisms
#' @aliases getAutomorphisms
#' @import GenomicRanges
#' @import IRanges
#' @importFrom methods new
#' @import S4Vectors
#' @export
#' @return This function returns an \code{\link{AutomorphismList-class}}
#' object as a list of \code{\link{Automorphism-class}} objects, which inherits
#' from \code{\link[GenomicRanges]{GRanges-class}} objects.
setMethod("getAutomorphisms",
    signature = "AutomorphismList",
    function(x) {
        nams <- names(x)
        gr <- x@SeqRanges
        x <- x@DataList
        if (length(gr) > 0) {
            x <- lapply(x, function(y) {
                if (inherits(y, c("DataFrame", "data.frame"))) {
                    mcols(gr) <- y
                    y <- as(gr, "Automorphism")
                }
                return(y)
            })
        } else {
            pos <- seq(1, length(x[[1]]), 1)
            gr <- GRanges(
                seqnames = 1,
                ranges = IRanges(start = pos, end = pos),
                strand = "+"
            )

            x <- lapply(x, function(y) {
                if (inherits(y, c("DataFrame", "data.frame"))) {
                    mcols(gr) <- y
                    y <- as(gr, "Automorphism")
                }
                return(y)
            })
        }
        x <- new("AutomorphismList",
            DataList = x,
            SeqRanges = gr
        )
        names(x) <- nams
        return(x)
    }
)


#' @rdname getAutomorphisms
#' @aliases getAutomorphisms
#' @import GenomicRanges
#' @export
#' @return An \code{\link{AutomorphismList-class}}
setMethod("getAutomorphisms",
    signature = "list",
    function(x) {
        x <- as.AutomorphismList(x)
        x <- getAutomorphisms(x)
        return(x)
    }
)


#' @rdname getAutomorphisms
#' @aliases getAutomorphisms
#' @import GenomicRanges
#' @export
#' @return An \code{\link{Automorphism-class}}
setMethod("getAutomorphisms",
    signature = "DataFrame_OR_data.frame",
    function(x) {
        as(x, "Automorphism")
    }
)
