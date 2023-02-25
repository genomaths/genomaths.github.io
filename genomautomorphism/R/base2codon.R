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

#' @rdname base2codon
#' @title Split a DNA sequence into codons
#' @description This function split a DNA sequence into a codon sequence.
#' @details It is expected that the provided DNA sequence is multiple of
#' 3, otherwise gaps are added to the end of the sequence.
#' @param x A character string, \code{\link[Biostrings]{DNAStringSet-class}}
#' or \code{\link[Biostrings]{DNAMultipleAlignment-class}} object carrying the
#' a DNA sequence.
#' @param ... Not in use.
#' @return If the argument of 'x' is character string, then a character vector
#' of codons will returned. If the argument of 'x' is
#' \code{\link[Biostrings]{DNAStringSet-class}} or
#' \code{\link[Biostrings]{DNAMultipleAlignment-class}} object, then a matrix
#' of codons is returned.
#' @author Robersy Sanchez <https://genomaths.com>. 01/15/2022
#' @examples
#'
#' ## Gaps are added at the sequence end.
#' seq <- c("ACCT")
#' base2codon(x = seq)
#'
#' ## This DNA sequence is multiple of 3
#' seq <- c("ACCTCA")
#' base2codon(x = seq)
#'
#' ## Load a DNAStringSet. A matrix of codons is returned
#' data(aln, package = "GenomAutomorphism")
#' base2codon(x = aln)
#'
#' @aliases base2codon
#' @export
setGeneric(
    "base2codon",
    function(x,
    ...) {
        standardGeneric("base2codon")
    }
)

#' @rdname base2codon
#' @aliases base2codon
#' @export
setMethod(
    "base2codon", signature(x = "character"),
    function(x) {
        nch <- nchar(x)
        if ((nch %% 3) != 0) {
            n <- (3 - nch %% 3)
            n <- paste(rep("-", n), collapse = "")
            x <- paste0(x, n)
            warning(
                "*** Base sequence of 'x' is not multiple of 3. ",
                "Gaps '-' have been added at the end of the sequence."
            )
        }
        n <- nchar(x)
        x <- substring(
            x,
            seq(1, n, 3),
            seq(3, n, 3)
        )
        return(x)
    }
)


#' @rdname base2codon
#' @aliases base2codon
#' @export
setMethod(
    "base2codon", signature(x = "DNAStringSet"),
    function(x) {
        x <- as.character(x)
        x <- lapply(x, base2codon)
        x <- do.call(rbind, x)
        return(t(x))
    }
)


#' @rdname base2codon
#' @aliases base2codon
#' @export
setMethod(
    "base2codon", signature(x = "DNAMultipleAlignment"),
    function(x) {
        x <- x@unmasked
        x <- base2codon(x)
        return(x)
    }
)
