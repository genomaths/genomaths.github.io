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

#' @rdname sortByChromAndStart
#'
#' @title Sorting \code{\link[GenomicRanges]{GRanges-class}} objects
#' @description Sorts a \code{\link[GenomicRanges]{GRanges-class}} objects
#' by seqname (chromosome), start, and position.
#' @details Objects that inherits from a
#' \code{\link[GenomicRanges]{GRanges-class}} can be sorted as well.
#' @param x GRanges object
#' @return \code{\link[GenomicRanges]{GRanges-class}} object or from the
#' original object class.
#' @examples
#' GR <- as(c("chr2:1-1", "chr1:1-1"), "GRanges")
#' GR <- sortByChromAndStart(GR)
#'
#' @importFrom BiocGenerics sort start end
#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqnames
#'
#' @aliases sortByChromAndStart
#' @export
sortByChromAndStart <- function(x) {
    seqlevels(x) <- sort(seqlevels(x))
    return(x[order(as.factor(seqnames(x)), start(x)), ])
}

#' @rdname sortByChromAndStart
#' @aliases sortByChromAndEnd
#' @export
sortByChromAndEnd <- function(x) {
    seqlevels(x) <- sort(seqlevels(x))
    return(x[order(as.factor(seqnames(x)), end(x)), ])
}
