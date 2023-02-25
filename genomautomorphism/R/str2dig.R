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


#' @aliases str2dig
#' @rdname str2dig
#' @title String to Digits
#' @description A simple function to transform a string of digits into a 
#' numeric vector.
#' @param x A character string or a list/ of character strings of 
#' numeric/digit symbols.
#' @param split The same as in \code{\link[base]{strsplit}}
#' @param num.cores,tasks Parameters for parallel computation using package
#' \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#' use, i.e. at most how many child processes will be run simultaneously (see
#' \code{\link[BiocParallel]{bplapply}} and the number of tasks per job (only
#' for Linux OS).
#' @param verbose If TRUE, prints the function log to stdout.
#' @param ... Further parameters for \code{\link[base]{strsplit}}.
#' @export
#' @returns A integer vector or a list of integer vectors.
#' @author Robersy Sanchez <https://genomaths.com>
#' @examples
#' ## A integer vector
#' str2dig("12231456247")
#'
#' ## A list of integer vectors
#' str2dig(list(num1 = "12231456247", num2 = "521436897"))


setGeneric("str2dig",
    function(
            x,
            split = "",
            ...) standardGeneric("str2dig"))


#' @aliases str2dig
#' @rdname str2dig
#' @export
setMethod("str2dig", signature(x = "character"),
    function(
        x, 
        split = "", 
        ...) {
        if (length(x) > 1) {
            x <- str2chr(x, split = "", ...)
            x <- lapply(x, as.integer)
        }
        else {
            x <- str2chr(x, split = "", ...)
            x <- as.integer(x)
        }
        return(x)
    }
)

#' @aliases str2dig
#' @rdname str2dig
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @export
setMethod("str2dig", signature(x = "list"),
    function(
        x,
        split = "",
        num.cores = 1L,
        tasks = 0L,
        verbose = FALSE,
        ...) {

        if (num.cores == 1) {
            x <- lapply(x, str2dig, split = split, ...)
        }
        else {
            ## ------------ Setting parallel computing ---------------- ##
            progressbar <- FALSE
            if (verbose) progressbar <- TRUE
            if (Sys.info()["sysname"] == "Linux") {
                bpparam <- MulticoreParam(
                    workers = num.cores,
                    tasks = tasks,
                    progressbar = progressbar)
            } else {
                bpparam <- SnowParam(
                    workers = num.cores,
                    type = "SOCK",
                    progressbar = progressbar)
            }
            ## --------------------------------------------------------- ##

            x <- bplapply(x, str2dig, split = split, ..., BPPARAM = bpparam)
        }
        return(x)
    }
)
