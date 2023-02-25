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

#' @rdname seqranges
#' @title Get DNA sequence Ranges and Coordinates representation on a given
#' Abelian Group
#' @description Extract the gene ranges and coordinates from a pairwise
#' alignment of codon/base sequences represented on a given Abelian group.
#' @param x An object from a \code{\link[Biostrings]{DNAStringSet}} or
#' \code{\link[Biostrings]{DNAMultipleAlignment}} class carrying the DNA
#' pairwise alignment of two sequences.
#' @param base_seq Logical. Whether to return the base or codon coordinates on
#' the selected Abelian group. If codon coordinates are requested, then the
#' number of the DNA bases in the given sequences must be multiple of 3.
#' @param granges Logical. Whether to return a
#' \code{\link[GenomicRanges]{GRanges-class}} object or a
#' \code{\link[S4Vectors]{DataFrame}}.
#' @param filepath A character vector containing the path to a file in
#' \emph{\strong{fasta}} format to be read. This argument must be given if
#' \emph{codon & base} arguments are not provided.
#' @param start,end,chr,strand Optional parameters required to build a
#' \code{\link[GenomicRanges]{GRanges-class}}. If not provided the default
#' values given for the function definition will be used.
#' @param ... Not in use.
#' @details This function provide an alternative way to get the codon
#' coordinate and the information on the codon sequence from a
#' \code{\link{CodonSeq}} class objects. The function can either take the
#' output from functions \code{\link{codon_coord}} or to operate directly on a
#' \code{\link[Biostrings]{DNAStringSet}} or to retrieve the a DNA sequence
#' alignment from a file.
#' @import S4Vectors
#' @importFrom methods new
#' @export
#' @seealso \code{\link{matrices}}, \code{\link{codon_coord}}, and
#' \code{\link{base_coord}}.
#' @author Robersy Sanchez <https://genomaths.com>
#' @references
#' \enumerate{
#'  \item Robersy Sanchez, Jesus Barreto (2021) Genomic Abelian Finite
#'   Groups.
#'  [doi:10.1101/2021.06.01.446543](https://doi.org/10.1101/2021.06.01.446543)
#'  \item M. V Jose, E.R. Morgado, R. Sanchez, T. Govezensky, The 24 possible
#'  algebraic representations of the standard genetic code in six or in three
#'  dimensions, Adv. Stud. Biol. 4 (2012) 119-152.[PDF](https://is.gd/na9eap).
#'  \item R. Sanchez. Symmetric Group of the Genetic-Code Cubes. Effect of the
#'  Genetic-Code Architecture on the Evolutionary Process MATCH Commun. Math.
#'  Comput. Chem. 79 (2018) 527-560.
#' }
#' @examples
#' ## Load a pairwise alignment
#' data(aln, package = "GenomAutomorphism")
#' aln
#'
#' ## A GRanges object carrying the aligned DNA sequence.
#' seqranges(
#'     x = aln,
#'     base_seq = TRUE,
#'     filepath = NULL,
#' )
#'
#' ## A GRanges object carrying the aligned codon sequence.
#' seqranges(
#'     x = aln,
#'     base_seq = FALSE,
#'     filepath = NULL,
#' )
#' @aliases seqranges
#' @return A \code{\link[GenomicRanges]{GRanges-class}}
setGeneric(
    "seqranges",
    function(x,
    ...) {
        standardGeneric("seqranges")
    }
)

#' @aliases seqranges
#' @rdname seqranges
#' @import GenomicRanges
#' @export
setMethod(
    "seqranges", signature(x = "CodonSeq"),
    function(x,
    granges = TRUE) {
        x <- x@SeqRanges
        if (granges) {
            x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
        }
        return(x)
    }
)


#' @aliases seqranges
#' @rdname seqranges
#' @export
setMethod(
    "seqranges", signature(x = "DNAStringSet_OR_NULL"),
    function(x,
    granges = TRUE,
    base_seq = TRUE,
    filepath = NULL,
    start = NA,
    end = NA,
    chr = 1L,
    strand = "+") {
        if (!is.null(filepath)) {
            x <- NULL
        }
        if (base_seq) {
            x <- base_coord(
                base = x,
                filepath = filepath,
                start = start,
                end = end,
                chr = chr,
                strand = strand
            )
        } else {
            x <- codon_coord(
                codon = x,
                filepath = filepath,
                group = "Z64",
                start = start,
                end = end,
                chr = chr,
                strand = strand
            )
        }
        x <- get_coord(x, output = "all")
        return(seqranges(x, granges = granges))
    }
)
