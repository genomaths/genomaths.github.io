## Copyright (C) 2022 Robersy Sanchez <https://genomaths.com/>
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

#' @rdname as.AutomorphismList
#' @aliases as.AutomorphismList
#' @title Methods for AutomorphismList-class Objects
#' @description Several methods are available to be applied on 
#' \code{\link{Automorphism-class}} and \code{\link{AutomorphismList-class}} 
#' objects.
#' @param x A \code{\link[S4Vectors]{DataFrame}} or a
#' \code{\link{automorphisms}} class object.
#' @param grs A \code{\link[GenomicRanges]{GRanges-class}} object.
#' @param ... Not in use yet.
#' @return The returned an AutomorphismList-class object.
#' @aliases as.AutomorphismList
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom methods setGeneric new
#' @export
#' @seealso \code{\link{automorphism_bycoef}}, \code{\link{automorphisms}}
#' @examples
#' ## Load a dataset
#' data("brca1_autm", package = "GenomAutomorphism")
#' 
#' ## Let's transforming into a list of Automorphisms-class objects
#' x1 <- as.list(brca1_autm[seq(2)])
#' 
#' ## Now, object 'x1' is transformed into a AutomorphismList-class object
#' as.AutomorphismList(x1)
#' 
#' ## Alternatively, let's transform the list 'x1' into a GRangesList-class 
#' ## object.
#' x1 <- GRangesList(x1)
#' 
#' ## Next, object 'x1' is transformed into a AutomorphismList-class object
#' as.AutomorphismList(x1)
setGeneric(
    "as.AutomorphismList",
    function(x,
            grs = GRanges(),
            ...) {
        standardGeneric("as.AutomorphismList")
    }
)


#' @rdname as.AutomorphismList
#' @aliases as.AutomorphismList
#' @import S4Vectors
#' @importFrom methods new
#' @export
setMethod(
    "as.AutomorphismList",
    signature(x = "GRangesList", grs = "GRanges_OR_NULL"),
    function(x,
            grs = GRanges(),
            ...) {
        if (length(grs) == 0) {
            grs <- x[[1]]
        }
        mcols(grs) <- NULL
        
        x <- lapply(x, function(y) {
            y <- as(y, "Automorphism")
            gr <- y
            mcols(gr) <- NULL
            if (any(gr != grs)) {
                stop("*** The ranges from the GRanges-class objects
                    must equals.")
            }
            return(mcols(y))
        })
        
        new("AutomorphismList",
            DataList = x,
            SeqRanges = grs
        )
    }
)

#' @rdname as.AutomorphismList
#' @aliases as.AutomorphismList
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom methods new
#' @export
setMethod(
    "as.AutomorphismList",
    signature(x = "list", grs = "GRanges_OR_NULL"),
    function(x,
            grs = GRanges(),
            ...) {
        if (length(grs) == 0) {
            if (inherits(x[[1]], "GRanges")) {
                grs <- x[[1]]
            } else {
                if (inherits(x[[1]], c("DataFrame", "data.frame"))) {
                    pos <- seq(1, nrow(x[[1]]), 1)
                    grs <- GRanges(
                        seqnames = 1,
                        ranges = IRanges(start = pos, end = pos),
                        strand = "+"
                    )
                } else {
                    stop(
                        "*** The argument of 'x' must be a list of ",
                        "objects from any of the classes: 'GRanges', ",
                        "'DataFrame', or 'data.frame'."
                    )
                }
            }
        }
        
        if (!is.null(mcols(grs))) {
            mcols(grs) <- NULL
        }
        
        if (all(slapply(x, function(y) inherits(y, "GRanges")))) {
            if (length(grs) == length(x)) {
                grs <- x
                mcols(grs) <- NULL
            }
            
            if (length(grs) == 0) {
                grs <- x[[1]]
                mcols(grs) <- NULL
            }
            
            x <- lapply(x, function(y) {
                y <- as(y, "Automorphism")
                return(mcols(y))
            })
            
            x <- new("AutomorphismList",
                    DataList = x,
                    SeqRanges = grs
            )
        }
        if (!inherits(x, "AutomorphismList")) {
            if (all(slapply(x, function(y) inherits(y, "DataFrame")))) {
                x <- new("AutomorphismList",
                        DataList = x,
                        SeqRanges = grs
                )
            }
        }
        return(x)
    }
)


## ======================== Show AutomorphismList ==================== #

#'@rdname AutomorphismList
#' @aliases show-AutomorphismList
#' @title Show method for \code{\link{AutomorphismList-class}} object
#' @param object An object from \code{\link{AutomorphismList-class}}.
#' @importFrom methods show
#' @import S4Vectors
#' @keywords internal
setMethod(
    "show",
    signature = "AutomorphismList",
    definition = function(object) {
        nams <- names(object@DataList)
        l <- length(nams)
        if (l > 10) {
            nams <- nams[c(seq(4), l - 2, l - 1, l)]
            nams[4] <- "..."
        }
        cat(class(object), " object of length: ",
            length(object@DataList), "\n",
            sep = ""
        )
        cat(paste0("names(", l, "):"), nams, "\n")
        cat("------- \n")
        gr <- object@SeqRanges
        if (length(gr) > 0 && inherits(object@DataList[[1]], "DataFrame")) {
            mcols(gr) <- object@DataList[[1]]
        } else {
            gr <- object@DataList[[1]]
        }
        
        print(as(gr, "Automorphism"))
        cat("...\n")
        cat("<", l - 1, " more ",
            class(object@DataList[[1]])[1], " element(s)>\n",
            sep = ""
        )
        cat("Two slots: 'DataList' & 'SeqRanges'\n")
        cat("------- \n")
        invisible(object)
    }
)


