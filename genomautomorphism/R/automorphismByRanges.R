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

#' @aliases automorphismByRanges
#' @rdname automorphismByRanges
#' @title Get the automorphisms by ranges.
#' @description Automorphisms estimated on a pairwise or a MSA alignment
#' can be grouped by ranges which inherits from
#' \code{\link[GenomicRanges]{GRanges-class}} or a
#' \code{\link[GenomicRanges]{GRanges-class}}.
#'
#' @param x An Automorphism-class object returned by function
#' \code{\link{automorphisms}}.
#' @param ... Not in use.
#' @return A  \code{\link[GenomicRanges]{GRanges-class}} or a
#' \code{\link[GenomicRanges]{GRangesList-class}}. Each
#' \code{\link[GenomicRanges]{GRanges-class}} object with a column
#' named *cube*, which carries the type of _cube_ automorphims.
#'
#' @export
#' @examples
#' ## Load dataset
#' data(autm, package = "GenomAutomorphism")
#'
#' automorphismByRanges(x = autm[c(1, 4)])
#'
setGeneric(
    "automorphismByRanges",
    function(x,
    ...) {
        standardGeneric("automorphismByRanges")
    }
)


#' @aliases automorphismByRanges
#' @rdname automorphismByRanges
#' @import GenomicRanges
#' @importFrom data.table data.table
#' @export
setMethod(
    "automorphismByRanges",
    signature(x = "Automorphism"),
    function(x) {
        end <- NULL
        i <- 1
        l <- length(x)

        if (!inherits(x, "GRanges")) {
            gr <- try(x@SeqRanges, silent = TRUE)
            if (!inherits(gr, "try-error")) {
                x <- getAutomorphisms(x)
                x <- x@DataList[[1]]
            }
        }

        idx <- vector(mode = "numeric", length = length(x))
        cube <- x$cube[1]
        for (k in seq_len(l)) {
            if (x$cube[k] != cube) {
                i <- i + 1
                cube <- x$cube[k]
            }
            idx[k] <- i
        }

        x$idx <- factor(idx)
        x <- data.table(data.frame(x))
        x <- x[, list(
            seqnames = unique(seqnames), start = min(start),
            end = max(end), strand = unique(strand),
            cube = unique(cube)
        ),
        by = idx
        ]
        x <- x[, c("seqnames", "start", "end", "strand", "cube")]
        x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
        x <- sortByChromAndStart(x)
        return(x)
    }
)


#' @aliases automorphismByRanges
#' @rdname automorphismByRanges
#' @param x An AutomorphismList-class object returned by function
#' \code{\link{automorphisms}}.
#' @param min.len Minimum length of a range to be reported.
#' @param num.cores,tasks Integers. Argument \emph{num.cores} denotes the
#' number of cores to use, i.e. at most how many child processes will be run
#' simultaneously (see \code{\link[BiocParallel]{bplapply}} function from
#' BiocParallel package). Argument \emph{tasks} denotes the number of tasks per
#' job. value must be a scalar integer >= 0L. In this documentation a job is
#' defined as a single call to a function, such as
#' \code{\link[BiocParallel]{bplapply}}. A task is the division of the \eqn{X}
#' argument into chunks. When tasks == 0 (default), \eqn{X} is divided as
#' evenly as possible over the number of workers (see
#' \code{\link[BiocParallel]{MulticoreParam}} from BiocParallel package).
#' @param verbose logic(1). If TRUE, enable progress bar.
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom parallel detectCores
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @importFrom data.table data.table
#' @export
setMethod(
    "automorphismByRanges", signature(x = "AutomorphismList"),
    function(x,
    min.len = 0L,
    num.cores = detectCores() - 1,
    tasks = 0L,
    verbose = TRUE) {
        gr <- try(x@SeqRanges, silent = TRUE)
        if (!inherits(gr, "try-error")) {
            x <- getAutomorphisms(x)
        }

        x <- x@DataList

        ## ---------------- Setting parallel computaton --------------- ##

        progressbar <- FALSE
        if (verbose) progressbar <- TRUE
        if (Sys.info()["sysname"] == "Linux") {
            bpparam <- MulticoreParam(
                workers = num.cores, tasks = tasks,
                progressbar = progressbar
            )
        } else {
            bpparam <- SnowParam(
                workers = num.cores, type = "SOCK",
                progressbar = progressbar
            )
        }

        ## -------------------------------------------------------------- ##

        if (length(gr) > 0) {
            x <- bplapply(x, function(x) {
                mcols(gr) <- x
                x <- automorphismByRanges(x)
                return(x)
            }, BPPARAM = bpparam)
        }

        idx <- which(slapply(x, function(x) {
            length(x) > min.len
        }))
        x <- x[idx]

        return(as(x, "GRangesList"))
    }
)
