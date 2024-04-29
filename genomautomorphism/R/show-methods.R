## Copyright (C) 2022-2024 Robersy Sanchez <https://genomaths.com/>
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

## ======================= Show - CodonSeq ==============================

#' @rdname show-methods
#' @aliases show-CodonSeq
#' @title Show method for 'CodonSeq' class object
#' @param object An object from 'CodonSeq'.
#' @importFrom methods show
#' @keywords internal
#' @export
setMethod(
    "show",
    signature = "CodonSeq",
    definition = function(object) {
        nams <- names(object@CoordList)
        cat(class(object), " object of length: ",
            length(object@CoordList), "\n",
            sep = ""
        )
        cat(paste0("names(", length(object@CoordList), "):"), nams, "\n")
        cat("------- \n")
        print(.showMatrix(object@CoordList[[1]]))
        cat("...\n")
        cat("<",
            length(object@CoordList) - 1,
            " more ", class(object@CoordList[[1]])[1], " element(s)>\n",
            sep = ""
        )
        cat("Two slots: 'CoordList' & 'SeqRanges'\n")
        cat("------- \n")
        invisible(object)
    }
)


## ======================= Show ListCodonMatrix =======================

#' @rdname show-methods
#' @aliases show-ListCodonMatrix
#' @title Show method for 'ListCodonMatrix' class object
#' @param object An object from 'ListCodonMatrix' class
#' @importFrom methods show
#' @keywords internal
#' @export
#' @return Print/show of a ListCodonMatrix-class object.
setMethod(
    "show",
    signature = "ListCodonMatrix",
    definition = function(object) {
        l <- length(object@DataList) 
        nms <- object@names
        cat(class(object), " object of length: ", l, "\n", sep = "")
        
        if (sum(nchar(nms)) > 50 && l > 1) 
            nms <- paste(nms[seq(6)])
        
        if (l > 1)
            cat("Seq.Alias:", nms, "...\n")
        else
            cat("Seq.Alias:", nms, "\n")
        
        cat("\n------- \n")
        
        cat("$", nms[1], "\n", sep = "")
        print(object[[1]])
        
        if (l > 2) {
            cat("\n\n")
            cat("$", nms[2], "\n", sep = "")
            print(object[[2]])
            cat("\n...\n")
            cat("\n<", l - 2, " more ", 
                class(object@DataList[[1]])[1], " element(s)>\n",
                sep = ""
            )
        }
        
        cat("Three slots: 'DataList', 'group', 'cube' & 'seq_alias'\n")
        cat("------- \n")
        invisible(object)
    }
)

## ========================= Show MatrixSeq ============================

#' @rdname MatrixSeq
#' @aliases show-MatrixSeq
#' @title Show method for 'MatrixSeq' class object
#' @param object An object from 'MatrixSeq' class
#' @importFrom methods show
#' @keywords internal
#' @export
#' @return Print/show of a MatrixSeq-class object.
setMethod(
    "show",
    signature = "MatrixSeq",
    function(object) {
        d <- dim(object@matrix)
        d2 <- d[2]
        if (d2 > 1)
            cat("MatrixSeq with", d[1], "rows (sequences) and", d2, 
            "columns (aminoacids/codons):\n")
        else {
            if (d[1] > 1)
                cat("MatrixSeq with", d[1], "rows (aminoacids/codons) and", 
                    d2, "column (sequence):\n")
        }
        cat("------- \n")
        if (d[1] <= 40 && d2 <= 20) {
            print(object@matrix)
        }
        else {
            x <- data.frame(object@matrix)
            cn <- colnames(x)
            rn <- rownames(x)
            
            if (d2 > 12) {
                x <- x[, c(seq(6), seq(d2 - 5, d2))]
                x[, 6] <- "..."
                cn <- colnames(x)
                cn[6] <- "..."
                colnames(x) <- cn
            }
            
            if (d[1] > 40) {
                x <- x[c(seq(11), seq(d[1] - 10, d[1])), ]
                if (is.numeric(x)) {
                    x <- data.frame(x)
                    rn <- rn[c(seq(11), seq(d[1] - 10, d[1]))] 
                    rn[11] <- "..."
                    rownames(x) <- rn
                    colnames(x) <- cn
                }
                else {
                    x[11, ] <- "..."
                }
            }
            print(x)
        }
        cat("------- \n")
        cat("Slots: 'seqs', 'matrix', 'names', 'aaindex',",
            "'phychem', 'accession")
        invisible(object)
    }
)


## ========================= Show MatrixList =====================

#' @rdname show-methods
#' @aliases show-MatrixList
#' @title Show method for 'MatrixList' class object
#' @param object An object from 'MatrixList' class
#' @importFrom methods show
#' @keywords internal
#' @export
#' @return Print/show of a MatrixList-class object.
setMethod(
    "show",
    signature = "MatrixList",
    definition = function(object) {
        nams <- names(object@matrices)
        cat(class(object), " object of length: ",
            length(object@matrices), "\n",
            sep = ""
        )
        cat(paste0("names(", length(object@matrices), "):"), nams, "\n")
        cat("------- \n")
        print(.showMatrix(object@matrices[[1]]))
        cat("...\n")
        cat("<",
            length(object@matrices) - 1,
            " more ", class(object@matrices[[1]])[1], " element(s)>\n",
            sep = ""
        )
        cat("Two slots: 'matrices' & 'names'\n")
        cat("------- \n")
        invisible(object)
    }
)

.showMatrix <- function(x) {
    d <- dim(x)
    if (!is.null(d)) {
        cat("Matrix with", d[1], "rows and", d[2], "columns:\n")
        if (d[1] > 10) {
            r <- c()
            for (k in c(seq(5), seq(d[1] - 5, d[1]))) {
                r <- rbind(r, x[k, ])
            }
            r[6, ] <- "..."
        }
        r <- data.frame(r)
        rown <- paste0(c(seq(5), seq(d[1] - 5, d[1])), ":")
        rown[6] <- "..."
        rownames(r) <- rown
    } else {
        l <- length(x)
        cat("Vector of length:", l, "\n")
        if (l > 10) {
            r <- x[c(seq(5), seq(l - 5, l))]
            r[6] <- "..."
            r <- data.frame(matrix(r, ncol = length(r)))
            nms <- paste0(c(seq(5), seq(l - 5, l)), ":")
            nms[6] <- r[6]
            colnames(r) <- nms
        }
        r <- x
    }
    return(r)
}

## ======================= Auxiliary functions ======================

#' @import S4Vectors
make_zero_col_DataFrame <- function(nrow = 0L) {
    stopifnot(isSingleNumber(nrow))
    if (!is.integer(nrow)) {
        nrow <- as.integer(nrow)
    }
    stopifnot(nrow >= 0L)
    new2("DFrame", nrows = nrow, check = FALSE)
}

#' @importFrom methods extends is
#' @import S4Vectors
new_SimpleList_from_list <- function(Class, x, type, ..., mcols) {
    if (!extends(Class, "SimpleList")) {
        stop("class ", Class, " must extend SimpleList")
    }
    if (!is.list(x)) {
        stop("'x' must be a list")
    }
    if (is.array(x)) {
        tmp_names <- names(x)
        dim(x) <- NULL
        names(x) <- tmp_names
    }
    class(x) <- "list"
    proto <- new(Class)
    if (missing(type)) {
        ans_elementType <- elementType(proto)
    } else {
        ans_elementType <- type
    }
    if (is(S4Vectors::mcols(proto, use.names = FALSE), "DataFrame")) {
        mcols <- make_zero_col_DataFrame(length(x))
    }
    if (!all(slapply(x, function(xi) extends(class(xi), ans_elementType)))) {
        stop(
            "all elements in 'x' must be ", ans_elementType,
            " objects."
        )
    }
    if (missing(mcols)) {
        return(new2(Class, listData = x, ..., check = FALSE))
    }
    new2(Class, listData = x, ..., elementMetadata = mcols, check = FALSE)
}



