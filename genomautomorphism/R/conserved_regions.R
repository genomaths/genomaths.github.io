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

#' @rdname conserved_regions
#' @aliases conserved_regions
#' @title Conserved and Non-conserved Regions from a MSA
#' @description Returns the Conserved or the Non-conserved Regions from a MSA.
#' @param x A \code{\link{Automorphism-class}}, a
#' \code{\link{AutomorphismList-class}},
#' a \code{\link{AutomorphismByCoef}} or a
#' \code{\link{AutomorphismByCoefList}} class object.
#' @param conserved Logical, Whether to return the \emph{conserved} or the
#' \emph{non-conserved regions}.
#' @param output A character string. Type of output.
#' @param ... Not in use.
#' @return A \code{\link{AutomorphismByCoef}} class object containing the
#' requested regions.
#' @import S4Vectors
#' @export
#' @examples
#' ## Load dataset
#' data("autm", package = "GenomAutomorphism")
#' conserved_regions(autm[1:3])
setGeneric(
    "conserved_regions",
    function(x,
    ...) {
        standardGeneric("conserved_regions")
    }
)


#' @rdname conserved_regions
#' @aliases conserved_regions
#' @import GenomicRanges
#' @export
setMethod("conserved_regions",
    signature = "Automorphism",
    function(x,
    conserved = TRUE,
    output = c("all_pairs", "unique_pairs", "unique")) {
        output <- match.arg(output)

        x <- automorphism_bycoef(x)
        x <- conserved_regions(
            x,
            conserved = conserved,
            output = output
        )
        return(x)
    }
)


#' @rdname conserved_regions
#' @aliases conserved_regions
#' @import GenomicRanges
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
#' @importFrom BiocParallel multicoreWorkers
#' @export
#' @examples
#' ## Load automorphism found COVID datatset
#' data("covid_autm", package = "GenomAutomorphism")
#' 
#' ## Conserved regions in the first 100 codons
#' conserv <- conserved_regions(covid_autm[1:100], output = "unique")
#' conserv
setMethod("conserved_regions",
    signature = "AutomorphismList",
    function(x,
    conserved = TRUE,
    output = c("all_pairs", "unique_pairs", "unique"),
    num.cores = multicoreWorkers(),
    tasks = 0L,
    verbose = FALSE) {
        output <- match.arg(output)

        x <- automorphism_bycoef(x,
            num.cores = num.cores,
            tasks = tasks,
            verbose = verbose
        )
        x <- unlist(x)
        x <- conserved_regions(
            x,
            conserved = conserved,
            output = output
        )
        return(x)
    }
)

#' @rdname conserved_regions
#' @aliases conserved_regions
#' @import GenomicRanges
#' @export
setMethod("conserved_regions",
    signature = "AutomorphismByCoef",
    function(x,
    conserved = TRUE,
    output = c("all_pairs", "unique_pairs", "unique")) {
        autm <- end <- NULL
        output <- match.arg(output)
        
        if (inherits(x$autm, "numeric")) {
            if (conserved) {
                x <- x[x$autm == 1]
            } else {
                x <- x[x$autm != 1]
            }
        }
        else {
            if (conserved) {
                x <- x[x$autm == "1,1,1"]
            } else {
                x <- x[x$autm != "1,1,1"]
            }
        }

        x <- sortByChromAndEnd(x)

        x <- switch(output,
            all_pairs = x,
            unique_pairs = unique(x),
            unique = {
                x <- unique(x)
                x <- data.table(data.frame(x))
                x <- x[, list(
                    seqnames = unique(seqnames),
                    end = min(end),
                    strand = unique(strand),
                    autm = unique(autm)
                ),
                by = c("start", "cube")
                ]
                x <- makeGRangesFromDataFrame(x,
                    keep.extra.columns = TRUE
                )
                x <- x[, c("autm", "cube")]
            }
        )
        return(x)
    }
)

#' @rdname conserved_regions
#' @aliases conserved_regions
#' @import GenomicRanges
#' @export
setMethod("conserved_regions",
    signature = "AutomorphismByCoefList",
    function(x,
    conserved = TRUE,
    output = c("all_pairs", "unique_pairs", "unique")) {
        output <- match.arg(output)

        x <- unlist(x)
        x <- conserved_regions(
            x,
            conserved = conserved,
            output = output
        )
        return(x)
    }
)
