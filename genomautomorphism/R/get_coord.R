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

#' @rdname get_coord
#' @title DNA base/codon sequence and coordinates represented on a given
#' Abelian group.
#' @description Given a string denoting a codon or base from the DNA (or RNA)
#' alphabet and a genetic-code Abelian group as given in reference (1), this
#' function returns an object from \code{\link{CodonGroup-class}} carrying the
#' DNA base/codon sequence and coordinates represented on the given Abelian
#' group.
#' @param x An object from a \code{\link{BaseGroup-class}},
#' \code{\link{CodonGroup-class}},
#' \code{\link[Biostrings]{DNAStringSet}} or
#' \code{\link[Biostrings]{DNAMultipleAlignment}} class carrying the DNA
#' pairwise alignment of two sequences. Objects from
#' \code{\link{BaseGroup-class}} and \code{\link{CodonGroup-class}} are
#' generated with functions: \code{\link{base_coord}} and
#' \code{\link{codon_coord}}, respectively.
#' @param base_seq Logical. Whether to return the base or codon coordinates on
#' the selected Abelian group. If codon coordinates are requested, then the
#' number of the DNA bases in the given sequences must be multiple of 3.
#' @param filepath A character vector containing the path to a file in
#' \emph{\strong{fasta}} format to be read. This argument must be given if
#' \emph{codon & base} arguments are not provided.
#' @param cube A character string denoting one of the 24 Genetic-code cubes,
#' as given in references (2 2 3).
#' @param group A character string denoting the group representation for the
#' given base or codon as shown in reference (1).
#' @param start,end,chr,strand Optional parameters required to build a
#' \code{\link[GenomicRanges]{GRanges-class}}. If not provided the default
#' values given for the function definition will be used.
#' @param output  See 'Value' section.
#' @param ... Not in use.
#' @details Symbols '-' and 'N' usually found in DNA sequence alignments to
#' denote gaps and missing/unknown bases are represented by the number: '-1'
#' on Z4 and '0' in Z5. In Z64 the symbol 'NA' will be returned for codons
#' including symbols '-' and 'N'.
#'
#' Although the \code{\link{CodonGroup-class}} object returned by
#' functions \code{\link{codon_coord}} and \code{\link{base_coord}} are useful
#' to store genomic information, the base and codon coordinates are not given
#' on them as numeric magnitudes. Function \code{\link{get_coord}} provides
#' the way to get the coordinates in a numeric object in object from and still
#' to preserve the base/codon sequence information.
#'
#' @import S4Vectors
#' @import Biostrings
#' @importFrom methods new
#' @return An object from \code{\link{CodonGroup-class}} class is returned
#' when \emph{output = 'all'}. This has two slots, the first one carrying a
#' list of matrices and the second one carrying the codon/base sequence
#' information. That is, if \emph{x} is an object from
#' \code{\link{CodonGroup-class}} class, then a list of matrices of codon
#' coordinate can be retrieved as x@CoordList and the information on the
#' codon sequence as x@SeqRanges.
#'
#' if \emph{output = 'matrix.list'}, then an object from
#' \code{\link{MatrixList}} class is returned.
#' @export
#' @aliases get_coord
#' @examples
#' ## Load a pairwise alignment
#' data("aln", package = "GenomAutomorphism")
#' aln
#'
#' ## DNA base representation in the Abelian group Z5
#' coord <- get_coord(
#'     x = aln,
#'     cube = "ACGT",
#'     group = "Z5"
#' )
#' coord ## A list of vectors
#'
#' ## Extract the coordinate list
#' coordList(coord)
#'
#' ## Extract the sequence list
#' seqRanges(coord)
#'
#' ## DNA codon representation in the Abelian group Z64
#' coord <- get_coord(
#'     x = aln,
#'     base_seq = FALSE,
#'     cube = "ACGT",
#'     group = "Z64"
#' )
#' coord
#'
#' ## Extract the coordinate list
#' coordList(coord)
#'
#' ## Extract the sequence list
#' seqRanges(coord)
#'
setGeneric("get_coord", function(x, ...) standardGeneric("get_coord"))

#' @aliases get_coord
#' @rdname get_coord
#' @import S4Vectors
#' @import GenomicRanges
#' @export
setMethod(
    "get_coord", signature(x = "BaseGroup_OR_CodonGroup"),
    function(x, output = c("all", "matrix.list")) {
        output <- match.arg(output)

        grp <- x@group
        gr <- x
        x <- mcols(x)
        nms <- colnames(x)
        idx <- grep("coord[[:digit:]]+", colnames(x))

        idx_seq <- grepl("seq[[:digit:]]+", colnames(x))
        if (any(idx_seq)) {
            gr <- gr[, which(idx_seq)]
            gr <- makeGRangesFromDataFrame(data.frame(gr),
                keep.extra.columns = TRUE
            )
        }

        if (grp == "Z5^3") {
            m <- lapply(idx, function(k) {
                m <- do.call(rbind, strsplit(x[, k], ","))
                m <- apply(m, 2, as.numeric)
            })
        } else {
            m <- lapply(idx, function(k) {
                m <- x[, k]
                as.numeric(m)
            })
        }
        names(m) <- nms[idx]
        res <- switch(output,
            all = new("CodonSeq", CoordList = m, SeqRanges = gr),
            matrix.list = new("MatrixList",
                matrices = m, names = nms
            )
        )
        return(res)
    }
)


#' @aliases get_coord
#' @rdname get_coord
#' @export
setMethod(
    "get_coord", signature(x = "DNAStringSet_OR_NULL"),
    function(x,
    output = c("all", "matrix.list"),
    base_seq = TRUE,
    filepath = NULL,
    cube = "ACGT",
    group = "Z4",
    start = NA,
    end = NA,
    chr = 1L,
    strand = "+") {
        if (!is.null(filepath)) {
            x <- NULL
        }
        if (base_seq) {
            x <- base_coord(
                base = x, filepath = filepath, cube = cube, group = group,
                start = start, end = end, chr = chr, strand = strand
            )
        } else {
            x <- codon_coord(
                codon = x, filepath = filepath, cube = cube, group = group,
                start = start, end = end, chr = chr, strand = strand
            )
        }
        x <- get_coord(x, output = output)
        return(x)
    }
)
