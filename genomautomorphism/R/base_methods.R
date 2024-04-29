## Copyright (C) 2021-2024 Robersy Sanchez <https://genomaths.com/>
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

## =========================== Docs ===========================

#' DNA Sequences Methods
#' 
#' @rdname base_methods
#' @aliases base_methods
#' @description 
#' 
#' ## Base coordinates on a given Abelian group representation
#' 
#' Given a string denoting a codon or base from the DNA (or RNA) alphabet,
#' function \strong{\emph{base_coord}} return the base coordinates in the
#' specify genetic-code Abelian group, as given in reference (1).
#' 
#' ## DNA sequences to \code{\link[GenomicRanges]{GRanges}} of bases. 
#' 
#' Function \strong{\emph{seq2granges}} transform an object from  
#' \code{\link[Biostrings]{DNAStringSet}}, 
#' \code{\link[Biostrings]{DNAMultipleAlignment-class}} or a character into
#' an object from [BaseSeq]. 
#' 
#' ## BaseSeq-class object to DNAStringSet-class object.
#' 
#' Function \strong{\emph{base_seq2string_set}} transforms an object from 
#' [BaseSeq] into an object from \code{\link[Biostrings]{DNAStringSet-class}}.
#' 
#' @details
#' 
#' ## Function 'base_coord'
#' 
#' Function \strong{\emph{base_coord}} is defined only for pairwise 
#' aligned sequences. Symbols "-" and "N" usually found in DNA sequence
#' alignments to denote gaps and missing/unknown bases are represented by the
#' number: '-1' on Z4 and '0' on Z5. In Z64 the symbol 'NA' will be returned
#' for codons including symbols "-" and "N".
#' 
#' ## Functions 'seq2granges' and 'base_seq2string_set'
#' 
#' For the sake of brevity the metacolumns from the object returned by 
#' function 'seq2granges' are named as 'S1', 'S2', 'S3', and so on. The
#' original DNA sequence alias are stored in the slot named 'seq_alias'.
#' (see examples).
#' 
#' @returns Depending on the function called, different object will be
#' returned:
#' 
#' ## Function 'base_coord'
#' 
#' This function returns a \code{\link{BaseGroup}} object
#' carrying the DNA sequence(s) and their respective coordinates in the
#' requested Abelian group of base representation (one-dimension, "Z4" or
#' "Z5"). Observe that to get coordinates in the set of of integer numbers
#' ("Z") is also possible but they are not defined to integrate a Abelian
#' group. These are just used for the further insertion the codon set in the
#' 3D space (R^3).
#'  
#' ## Function 'seq2granges'
#' 
#' This function returns a [BaseGroup] object carrying the DNA sequence(s), one
#' base per ranges. A [BaseGroup] class object inherits from  
#' \code{\link[GenomicRanges]{GRanges-class}}.
#' 
#' ## Function 'base_seq2string_set'
#' 
#' This function returns a \code{\link[Biostrings]{DNAStringSet-class}}.
#' 

## ===================== base_coord ======================

#' @param base An object from a \code{\link[Biostrings]{DNAStringSet}} or
#' \code{\link[Biostrings]{DNAMultipleAlignment}} class carrying the DNA
#' pairwise alignment of two sequences.
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
#' 
#' @seealso [Symmetric Group of the Genetic-Code Cubes.](
#' https://github.com/genomaths/GenomeAlgebra_SymmetricGroup)
#' @import S4Vectors
#' @import Biostrings
#' @importFrom methods new
#' @export
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
#' @seealso \code{\link{codon_coord}} and \code{\link{base2int}}.
#' @examples
#' ## Example 1. Let's get the base coordinates for codons "ACG"
#' ## and "TGC":
#' x0 <- c("ACG", "TGC")
#' x1 <- DNAStringSet(x0)
#' x1
#' 
#' ## Get the base coordinates on cube = "ACGT" on the Abelian group = "Z4"
#' base_coord(x1, cube = "ACGT", group = "Z4")
#' 
#' ## Example 2. Load a pairwise alignment
#' data("aln", package = "GenomAutomorphism")
#' aln
#'
#' ## DNA base representation in the Abelian group Z4
#' bs_cor <- base_coord(
#'     base = aln,
#'     cube = "ACGT"
#' )
#' bs_cor
#'
#' ## Example 3. DNA base representation in the Abelian group Z5
#' bs_cor <- base_coord(
#'     base = aln,
#'     cube = "ACGT",
#'     group = "Z5"
#' )
#' bs_cor
#' 
#' ## Example 4. Load a multiple sequence alignment (MSA) of primate BRCA1 DNA  
#' ## repair genes 
#' data("brca1_aln2", package = "GenomAutomorphism")
#' brca1_aln2
#' 
#' ## Get BaseSeq-class object
#' gr <- seq2granges(brca1_aln2)
#' gr
#' 
#' ## Transform the BaseSeq-class object into a DNAStringSet-class object
#' str_set <- base_seq2string_set(gr)
#' str_set
#' 
#' ## Recovering the original MSA
#' DNAMultipleAlignment(as.character(str_set))
#' 
#' ## Example 5. 
#' base_matrix(base = aln, cube = "CGTA", group = "Z5")
#' 
#' @aliases base_coord
#' @export
#' @return A BaseGroup-class object.
setGeneric(
    "base_coord",
    function(base = NULL,
    filepath = NULL,
    cube = "ACGT",
    group = "Z4",
    ...) {
        standardGeneric("base_coord")
    }
)


#' @aliases base_coord
#' @rdname base_methods
#' @import GenomicRanges
#' @importFrom methods new
#' @import Biostrings
#' @export
setMethod(
    "base_coord", signature(base = "DNAStringSet_OR_NULL"),
    function(
        base = NULL,
        filepath = NULL,
        cube = c(
            "ACGT", "AGCT", "TCGA", "TGCA", "CATG",
            "GTAC", "CTAG", "GATC", "ACTG", "ATCG",
            "GTCA", "GCTA", "CAGT", "TAGC", "TGAC",
            "CGAT", "AGTC", "ATGC", "CGTA", "CTGA",
            "GACT", "GCAT", "TACG", "TCAG"),
        group = c("Z4", "Z5"),
        start = NA,
        end = NA,
        chr = 1L,
        strand = "+") {

        cube <- toupper(cube)
        group <- toupper(group)
        cube <- match.arg(cube)
        group <- match.arg(group)

        if (is.null(base) && is.null(filepath)) {
            stop(
                "*** Arguments 'base' & 'filepath' cannot be",
                " simultaneously NULL."
            )
        }

        if (!is.null(filepath) && is.character(filepath)) {
            base <- readDNAMultipleAlignment(filepath = filepath)
        }

        if (inherits(base, "DNAMultipleAlignment")) {
            base <- base@unmasked
        }

        if (length(base) > 2)
            stop("*** Function 'base_coord' is defined only for pairwise
                aligned sequences.")
            
        
        len <- min(width(base))

        if (!is.na(start) || !is.na(end)) {
            if (!is.na(start) && start > len) {
                stop(
                    "*** The 'start' argument is greater than",
                    " the 'base' length"
                )
            }
            if (!is.na(end) && end > len) {
                stop(
                    "*** The 'end' argument is greater than",
                    "   the 'base' length"
                )
            }
            base <- DNAStringSet(
                base,
                start = start,
                end = end
            )
        }

        if (length(base) > 1) {
            if (length(unique(width(base))) > 1) {
                stop("*** Function 'base_coord' is defined only for pairwise
                aligned sequences:\n", "'base' strings are not equal-width.")
            }
                
            base <- t(as.matrix(base))
            colnames(base) <- NULL
            base <- data.frame(base)
        } else {
            base <- as.character(base)
            base <- strsplit(base, "")[[1]]
            base <- data.frame(base)
        }
        seq <- base
        base <- base2int(base = base, cube = cube, group = group)
        if (is.na(start)) {
            start <- 1L
        }
        if (is.na(end)) {
            end <- len
        }

        pos <- seq(start, end, 1L)
        if (!is.null(dim(base))) {
            colnames(base) <- paste0("coord", seq_len(ncol(base)))
            colnames(seq) <- paste0("seq", seq_len(ncol(base)))
        }
        base <- data.frame(
            seqnames = chr,
            start = pos,
            end = pos,
            strand = strand,
            seq,
            base
        )

        base <- makeGRangesFromDataFrame(base, keep.extra.columns = TRUE)
        base <- new(
            "BaseGroup",
            seqnames = seqnames(base),
            ranges = ranges(base),
            strand = strand(base),
            elementMetadata = base@elementMetadata,
            seqinfo = base@seqinfo,
            colnames = colnames(base@elementMetadata),
            group = group,
            cube = cube
        )
        return(base)
    }
)

## ===================== seq2granges ======================

#' @rdname base_methods
#' @aliases seq2granges
#' @title DNA sequences to GRanges of bases and the reverse.
#' @seealso [Symmetric Group of the Genetic-Code Cubes.](
#' https://github.com/genomaths/GenomeAlgebra_SymmetricGroup)
#' @import S4Vectors
#' @import Biostrings
#' @importFrom methods new
#' @export
#' @author Robersy Sanchez <https://genomaths.com>
#' @seealso [base_coord] and [codon_coord].
#' @export
setGeneric(
    "seq2granges",
    function(
        base = NULL,
        filepath = NULL,
        ...) {
        standardGeneric("seq2granges")
    }
)


#' @aliases seq2granges
#' @rdname base_methods
#' @param seq_alias DNA sequence alias/ID and description.
#' @param ... Not in use yet.
#' @import GenomicRanges
#' @importFrom methods new
#' @import Biostrings
#' @export
setMethod(
    "seq2granges", signature(base = "DNAStringSet_OR_NULL"),
    function(
        base = NULL,
        filepath = NULL,
        start = NA,
        end = NA,
        chr = 1L,
        strand = "+",
        seq_alias = NULL,
        ...) {

    if (is.null(base) && is.null(filepath)) {
        stop(
            "*** Arguments 'base' & 'filepath' cannot be",
            " simultaneously NULL."
        )
    }
    
    if (!is.null(filepath) && is.character(filepath)) {
        base <- readDNAMultipleAlignment(filepath = filepath)
    }
    
    if (inherits(base, "DNAMultipleAlignment")) 
        base <- base@unmasked
    
    len <- min(width(base))
    
    if (!is.na(start) || !is.na(end)) {
        if (!is.na(start) && start > len) {
            stop(
                "*** The 'start' argument is greater than",
                " the 'base' length"
            )
        }
        if (!is.na(end) && end > len) {
            stop(
                "*** The 'end' argument is greater than",
                "   the 'base' length"
            )
        }
        base <- DNAStringSet(
            base,
            start = start,
            end = end
        )
    }
    
    if (inherits(base, "DNAStringSet") && is.null(seq_alias))
        seq_alias <- names(base)
    
    if (length(base) > 1) {
        base <- t(as.matrix(base))
        colnames(base) <- NULL
        base <- data.frame(base)
    } else {
        base <- as.character(base)
        base <- strsplit(base, "")[[1]]
        base <- data.frame(base)
    }
    seq <- base
    if (is.na(start)) {
        start <- 1L
    }
    if (is.na(end)) {
        end <- len
    }
    
    pos <- seq(start, end, 1L)
    if (!is.null(dim(base))) 
        colnames(base) <- paste0("S", seq_len(ncol(base)))
    
    base <- data.frame(
        seqnames = chr,
        start = pos,
        end = pos,
        strand = strand,
        base
    )
    
    base <- makeGRangesFromDataFrame(base, keep.extra.columns = TRUE)
    
    base <- new(
        "BaseSeq",
        seqnames = seqnames(base),
        ranges = ranges(base),
        strand = strand(base),
        elementMetadata = base@elementMetadata,
        seqinfo = base@seqinfo,
        seq_alias = seq_alias
    )
    
    return(base)
    }
)

## ===================== base_seq2string_set ======================


#' @aliases base_seq2string_set
#' @rdname base_methods
#' @import GenomicRanges
#' @importFrom methods new
#' @import Biostrings
#' @param x A 'BaseSeq' class object. 
#' @export
#' @examples
#' ## Example 5. 
#' 
setGeneric(
    "base_seq2string_set",
    function(
        x,
        ...) {
        standardGeneric("base_seq2string_set")
    }
)


#' @aliases base_seq2string_set
#' @rdname base_methods
#' @import GenomicRanges
#' @import Biostrings
#' @export
setMethod(
    "base_seq2string_set", signature(x = "BaseSeq"),
    function(x) {
        seq_alias <- x@seq_alias
        x <- mcols(x)
        x <- apply(x, 2, paste0, collapse = "")
        x <- DNAStringSet(x)
        names(x) <- seq_alias
        return(x)
    }
)


## ===================== base_matrix ======================

#' @aliases base_matrix
#' @rdname base_methods
#' @import GenomicRanges
#' @importFrom methods new
#' @import Biostrings
#' @export
setGeneric(
    "base_matrix",
    function(
        base,
        ...) {
        standardGeneric("base_matrix")
    }
)


#' @aliases base_matrix
#' @rdname base_methods
#' @import GenomicRanges
#' @import Biostrings
#' @export
setMethod(
    "base_matrix", signature(base = "DNAStringSet_OR_NULL"),
    function(
        base,
        cube = c(
            "ACGT", "AGCT", "TCGA", "TGCA", "CATG",
            "GTAC", "CTAG", "GATC", "ACTG", "ATCG",
            "GTCA", "GCTA", "CAGT", "TAGC", "TGAC",
            "CGAT", "AGTC", "ATGC", "CGTA", "CTGA",
            "GACT", "GCAT", "TACG", "TCAG"),
        group = c("Z4", "Z5"),
        seq_alias = NULL) {
        
        cube <- toupper(cube)
        group <- toupper(group)
        cube <- match.arg(cube)
        group <- match.arg(group)
        
        if (is.null(seq_alias) && inherits(base, "DNAStringSet"))
            seq_alias <- names(DNAStringSet)
        
        base <- seq2granges(base)
        if (is.null(seq_alias)) 
            seq_alias <- base@seq_alias
        
        mcols(base) <- base2int(base = data.frame(mcols(base)), 
                                cube = cube, group = group)
        
        base <- new(
            "BaseSeqMatrix",
            seqnames = seqnames(base),
            ranges = ranges(base),
            strand = strand(base),
            elementMetadata = base@elementMetadata,
            seqinfo = base@seqinfo,
            seq_alias = seq_alias,
            group = group,
            cube = cube
        )
        
        return(base)
    }
)



### -------------------------- Auxiliary functions ------------------------

#' Check URLs
#' @details Internal use only.
#' @keywords internal
#' @return Logical values
is.url <- function(x) {
    heads <- c(
        "//", "http://", "https://", "ftp://", "ftps://",
        "file://"
    )
    any(slapply(heads, function(p) grepl(pattern = p, x)))
}

