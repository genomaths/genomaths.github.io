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

## ========================== GRanges_OR_NULL =====================

#' @rdname GRanges_OR_NULL
#' @title A definition for the union of 'GRanges' and 'NULL' class.
#' @import GenomicRanges
#' @keywords internal
#' @export
setClassUnion("GRanges_OR_NULL", c("GRanges", "NULL", "missing"))


## ========================== BaseSeq =============================

setClassUnion("character_OR_NULL", c("character", "NULL", "missing"))

#' @aliases BaseSeq
#' @rdname BaseSeq
#' @title A class definition to store DNA base sequence.
#' @import S4Vectors
#' @import GenomicRanges
#' @seealso \code{\link{automorphisms}}
#' @keywords internal
#' @export
#' @return Given the slot values define a BaseSeq-class.
setClass("BaseSeq",
    slots = c(
        seqnames = "Rle",
        ranges = "IRanges_OR_IPos",
        strand = "Rle",
        elementMetadata = "DataFrame",
        seqinfo = "Seqinfo",
        seq_alias = "character_OR_NULL"
    ),
    contains = "GRanges"
)


## ========================== BaseSeqMatrix =========================


#' @aliases BaseSeqMatrix
#' @rdname BaseSeqMatrix
#' @title A class definition to Store DNA base sequence coordinates in a given
#' Genetic Code Cube.
#' @import S4Vectors
#' @import GenomicRanges
#' @seealso \code{\link{automorphisms}}
#' @keywords internal
#' @export
#' @return Given the slot values define a BaseSeq-class.
setClass("BaseSeqMatrix",
    slots = c(
        seqnames = "Rle",
        ranges = "IRanges_OR_IPos",
        strand = "Rle",
        elementMetadata = "DataFrame",
        seqinfo = "Seqinfo",
        seq_alias = "character_OR_NULL",
        group = "character",
        cube = "character"
    ),
    contains = "GRanges"
)


## ========================== BaseGroup =============================

#' @aliases BaseGroup
#' @rdname BaseGroup
#' @title A class definition to store codon automorphisms in given in the
#' Abelian group representation.
#' @import S4Vectors
#' @import GenomicRanges
#' @seealso \code{\link{automorphisms}}
#' @keywords internal
#' @export
#' @return Given the slot values define a BaseGroup-class.
setClass("BaseGroup",
    slots = c(
        seqnames = "Rle",
        ranges = "IRanges_OR_IPos",
        strand = "Rle",
        elementMetadata = "DataFrame",
        seqinfo = "Seqinfo",
        colnames = "character",
        group = "character",
        cube = "character"
    ),
    contains = "GRanges"
)

# ====================  Validity BaseGroup ======================== #

#' @rdname valid.BaseGroup
#' @aliases valid.BaseGroup.elem
#' @title Valid BaseGroup mcols
#' @param x A 'BaseGroup' object
#' @keywords internal
#' @return If valid return NULL
valid.BaseGroup.elem <- function(x) {
    m1 <- paste0(
        "*** This is not a valid  BaseGroup-class object.",
        " Columns from the matacolumn have the wrong names"
    )
    m2 <- paste0(
        "*** This is not a valid  BaseGroup-class object.",
        " seq1 or seq2 columns is not a base sequence"
    )
    m3 <- paste0(
        "*** Argument 'x' is not a BaseGroup-class object.",
        " The slot 'group' is not present or wrong naming."
    )
    m4 <- paste0(
        "*** Argument 'x' is not a BaseGroup-class object.",
        " The slot 'cube' is not present or wrong naming."
    )
    r1 <- r2 <- r3 <- r4 <- FALSE
    if (length(x) > 0) {
        coln <- x@colnames
        if (any(!is.element(
            coln,
            c("seq1", "seq2", "coord1", "coord2")
        ))) {
            r1 <- TRUE
        }
        if (unique(nchar(x$seq1)) != 1 || unique(nchar(x$seq2)) != 1) {
            r2 <- TRUE
        }

        group <- try(x@group, silent = TRUE)
        elem <- is.element(group, c("Z4", "Z5", "Z4^3", "Z5^3"))
        if (inherits(group, "try-error") || !elem) {
            r3 <- TRUE
        }

        cube <- try(x@cube, silent = TRUE)
        celem <- is.element(cube, c(
            "ACGT", "AGCT", "TCGA", "TGCA", "CATG",
            "GTAC", "CTAG", "GATC", "ACTG", "ATCG",
            "GTCA", "GCTA", "CAGT", "TAGC", "TGAC",
            "CGAT", "AGTC", "ATGC", "CGTA", "CTGA",
            "GACT", "GCAT", "TACG", "TCAG"
        ))
        if (inherits(celem, "try-error") || !celem) {
            r4 <- TRUE
        }
    }
    if (any(c(r1, r2, r3, r4))) {
        res <- c(m1, m2, m3, m4)[c(r1, r2, r3, r4)]
        return(res[1])
    }
    NULL
}

#' @rdname valid.BaseGroup
#' @aliases valid.GRanges
#' @title Valid 'BaseGroup' inheritance from 'GRanges' class
#' @param x A 'BaseGroup object'
#' @keywords internal
#' @return If valid return NULL
valid.GRanges <- function(x) {
    if (length(x) > 0) {
        if (!inherits(x, "GRanges")) {
            return("*** This is not a valid  BaseGroup-class object.")
        }
    }
    NULL
}

#' @rdname valid.BaseGroup
#' @aliases valid.BaseGroup
#' @title Valid BaseGroup
#' @param x A 'BaseGroup object'
#' @keywords internal
#' @return If valid return NULL
valid.BaseGroup <- function(x) {
    c(valid.GRanges(x), valid.BaseGroup.elem(x))
}

setValidity2("BaseGroup", valid.BaseGroup)


## ========================== CodonGroup =============================

#' @aliases CodonGroup
#' @rdname CodonGroup
#' @title A class definition to store codon automorphisms in given in the
#' Abelian group representation.
#' @import S4Vectors
#' @import GenomicRanges
#' @seealso \code{\link{automorphisms}}
#' @keywords internal
#' @export
#' @return Given the slot values define a CodonGroup-class.
setClass("CodonGroup",
    slots = c(
        seqnames = "Rle",
        ranges = "IRanges_OR_IPos",
        strand = "Rle",
        elementMetadata = "DataFrame",
        seqinfo = "Seqinfo",
        colnames = "character",
        group = "character",
        cube = "character"
    ),
    contains = "GRanges"
)

#' @rdname BaseGroup_OR_CodonGroup
#' @aliases BaseGroup_OR_CodonGroup
#' @title A definition for the union of classes 'BaseGroup' and  'CodonGroup'
#' @keywords internal
#' @seealso \code{\link{BaseGroup}} and \code{\link{CodonGroup}}.
#' @export
setClassUnion("BaseGroup_OR_CodonGroup", c("BaseGroup", "CodonGroup"))


# ====================  Validity CodonGroup ======================== #


#' @rdname valid.CodonGroup
#' @aliases valid.CodonGroup.mcols
#' @title Valid CodonGroup mcols
#' @param x A 'CodonGroup' object
#' @keywords internal
#' @return If valid return NULL
valid.CodonGroup.mcols <- function(x) {
    if (length(x) > 0) {
        coln <- x@colnames
        m1 <- paste0(
            "*** This is not a valid  CodonGroup-class object.",
            "Columns from the matacolumn have the wrong names"
        )
        m2 <- paste0(
            "*** This is not a valid  BaseGroup-class object.",
            "seq1 or seq2 columns is not a base-triplet sequence"
        )
        m3 <- paste0(
            "*** This is not a CodonGroup-class object.",
            " The slot 'cube' is not present or wrong naming."
        )
        m4 <- paste0(
            "*** Argument 'x' is not a CodonGroup-class object.",
            "The slot 'group' is not present or wrong naming."
        )
        r1 <- r2 <- r3 <- r4 <- FALSE

        if (any(!is.element(
            coln,
            c("seq1", "seq2", "coord1", "coord2")
        ))) {
            r1 <- TRUE
        }
        if (unique(nchar(x$seq1)) != 3 || unique(nchar(x$seq2)) != 3) {
            r2 <- TRUE
        }

        cube <- try(x@cube, silent = TRUE)
        celem <- is.element(cube, c(
            "ACGT", "AGCT", "TCGA", "TGCA", "CATG",
            "GTAC", "CTAG", "GATC", "ACTG", "ATCG",
            "GTCA", "GCTA", "CAGT", "TAGC", "TGAC",
            "CGAT", "AGTC", "ATGC", "CGTA", "CTGA",
            "GACT", "GCAT", "TACG", "TCAG"
        ))
        if (inherits(celem, "try-error") || !celem) {
            r3 <- TRUE
        }

        group <- try(x@group, silent = TRUE)
        elem <- is.element(group, c("Z4", "Z5", "Z4^3", "
                                    Z5^3", "Z64", "Z125"))
        if (inherits(group, "try-error") || !elem) {
            r4 <- TRUE
        }
    }
    if (any(c(r1, r2, r3, r4))) {
        res <- c(m1, m2, m3, m4)[c(r1, r2, r3, r4)]
        return(res[1])
    }
    NULL
}

#' @rdname valid.CodonGroup
#' @aliases valid.CodonGroup
#' @title Valid CodonGroup
#' @param x A 'CodonGroup object'
#' @keywords internal
#' @return If valid return NULL
valid.CodonGroup <- function(x) {
    c(valid.GRanges(x), valid.CodonGroup.mcols(x))
}

setValidity2("CodonGroup", valid.CodonGroup)


## ========================== CodonSeq =============================

#' @rdname CodonSeq
#' @aliases CodonSeq
#' @title A class definition to store codon coordinates given in the Abelian
#' group and the codon sequence.
#' @description An objects from 'CodonSeq' or 'MatrixList' class is returned 
#' by function \code{\link{get_coord}}. This object will store the coordinate
#' of each sequence in a list of 3D-vectors or a list of vectors located in 
#' the slot named 'CoordList'. The original codon sequence (if provided) will
#' be stored in the slot named 'SeqRanges'.
#' @importFrom methods validObject
#' @import GenomicRanges
#' @keywords internal
#' @export
#' @return Given the slot values define a CodonSeq-class.
setClass("CodonSeq",
    slots = c(
        CoordList = "list",
        SeqRanges = "GenomicRanges_OR_missing"
    )
)

#' @aliases coordList
#' @rdname CodonSeq
#' @title Method to extract 'CoordList' slot from a
#' \code{\link{CodonSeq-class}}
#' @param x An object from \code{\link{CodonSeq-class}}.
setGeneric(
    "coordList",
    function(x) standardGeneric("coordList")
)

#' @aliases coordList
#' @rdname CodonSeq
#' @export
#' @examples 
#' ## Load a DNA sequence alignment 
#' data("aln", package = "GenomAutomorphism")
#' 
#' ## Get base coordinates on 'Z5'
#' coord <- get_coord(
#'     x = aln,
#'     cube = "ACGT",
#'     group = "Z5"
#' )
#' coordList(coord)
setMethod(
    "coordList", signature(x = "CodonSeq"),
    function(x) x@CoordList
)

#' @aliases seqRanges
#' @rdname CodonSeq
#' @title Method to extract 'SeqRanges' slot from a
#' \code{\link{CodonSeq-class}}
#' @param x An object from \code{\link{CodonSeq-class}}.
setGeneric(
    "seqRanges",
    function(x) standardGeneric("seqRanges")
)

#' @aliases seqRanges
#' @rdname CodonSeq
#' @export
#' @examples 
#' ## Load a DNA sequence alignment
#' data("aln", package = "GenomAutomorphism")
#' 
#' ## Get base coordinates on 'Z5'
#' coord <- get_coord(
#'     x = aln,
#'     cube = "ACGT",
#'     group = "Z5"
#' )
#' 
#' seqRanges(coord)
setMethod(
    "seqRanges", signature(x = "CodonSeq"),
    function(x) x@SeqRanges
)

## ======================= CodonMatrix class =========================


#' @aliases CodonMatrix
#' @rdname CodonMatrix
#' @title A Convenient Class to Store a Codon Coordinate in given 
#' Genetic Code cube.
#' @description
#' A CodonMatrix is the object created by function [codon_matrix]
#' 
#' @keywords internal
#' @export
#' @examples
#' ## CodonMatrix is generated by function 'codon_matrix' (inside a
#' ## ListCodonMatrix-class object)
#' ## Let's create DNAStringSet-class object
#' base <- DNAStringSet(x = c(S1 = 'ACGTGATCAAGT'))
#' 
#' x1 <- codon_matrix(base)
#' x1
#' 
#' ## Extract the first element
#' x1[1]
#' x1$codon.1
#' x1[[1]]
setClass("CodonMatrix",
    slots = c(
        seqnames = "Rle",
        ranges = "IRanges_OR_IPos",
        strand = "Rle",
        elementMetadata = "DataFrame",
        seqinfo = "Seqinfo",
        group = "character",
        cube = "character",
        seq_alias = "character_OR_NULL"
    ),
    contains = "GRanges"
)

#' Constructor of CodonMatrix-class
#' @aliases CodonMatrix
#' @rdname CodonMatrix
#' @param object A \code{\link[GenomicRanges]{GRanges-class}} object.
#' @param group The name of the base group, 'Z4' or 'Z5'.
#' @param cube The name of the genetic-code cube applied to get the 
#' codon coordinates.
#' @param seq_alias The 'alias' given to the codon sequence.
#' @seealso [base_coord] and [codon_coord].
#' @returns A 'CodonMatrix' class object
#' @export
CodonMatrix <- function(object, group, cube, seq_alias = NULL) {
    new("CodonMatrix",
        seqnames = object@seqnames,
        ranges = object@ranges,
        strand = object@strand,
        elementMetadata = object@elementMetadata,
        seqinfo = object@seqinfo,
        group = group,
        cube = cube,
        seq_alias = seq_alias)
}

## ======================= ListCodonMatrix class =========================

#' @aliases ListCodonMatrix
#' @rdname ListCodonMatrix
#' @title A Convenient Class to Store Codon Coordinates in given 
#' Genetic Code cube.
#' @description
#'  ListCodonMatrix-class objects are generated by function [codon_matrix].
#' @keywords internal
#' @export
#' @examples
#' ## ListCodonMatrix-class objects are generated by function 'codon_matrix'.
#' ## Let's create DNAStringSet-class object
#' base <- DNAStringSet(x = c( seq1 ='ACGTGATCAAGT',
#'                             seq2 = 'GTGTGATCCAGT',
#'                             seq3 = 'TCCTGATCAGGT'))
#' 
#' x1 <- codon_matrix(base)
#' x1
#' 
#' ## Extract the first element
#' x1[1]
#' x1$codon.1
#' x1[[1]]
#'   
setClass("ListCodonMatrix",
    slots = c(
        DataList = "list",
        group = "character",
        cube = "character",
        seq_alias = "character_OR_NULL",
        names = "character_OR_NULL"
    )
)


#' Constructor of ListCodonMatrix-class
#' @aliases ListCodonMatrix
#' @rdname ListCodonMatrix
#' @param object A list of CodonMatrix-class objects
#' @export
ListCodonMatrix <- function(
    object, cube, group, seq_alias = NULL, names = NULL) {
    
    if (is.null(names))
        names <- paste0("codon.", seq_along(object))
    
    if (is.null(seq_alias))
        names <- paste0("seq", seq_along(object))
    
    new("ListCodonMatrix", 
        DataList = object, 
        cube = cube, 
        group = group, 
        seq_alias = seq_alias,
        names = names)
}


## ==================== Validity ListCodonMatrix ==================== #
#
#' @aliases valid.ListCodonMatrix
#' @rdname ListCodonMatrix
#' @title Valid ListCodonMatrix-class object 
#' @param x A 'ListCodonMatrix-class' object
#' @import S4Vectors
#' @keywords internal

valid.ListCodonMatrix <- function(x) {
    
    group <- x@group
    cube <- x@cube
    alias <- x@seq_alias
    
    r1 <- c("*** This is not a valid ListCodonMatrix class object.\n",
            "Every element from the list must be a CodonMatrix object.")
    r2 <- "*** The provided base-group must be 'Z4' or 'Z5'."
    r3 <- "*** The provided genetic-code cube name is not valid."
    
    t1 <- t2 <- t3 <- t4 <- FALSE
    
    t1 <- !all(slapply(x@DataList, inherits, what = "CodonMatrix"))
    t2 <- !is.element(group, c("Z5", "Z4"))
    t3 <- !is.element(cube, c(
        "ACGT", "AGCT", "TCGA", "TGCA", "CATG",
        "GTAC", "CTAG", "GATC", "ACTG", "ATCG",
        "GTCA", "GCTA", "CAGT", "TAGC", "TGAC",
        "CGAT", "AGTC", "ATGC", "CGTA", "CTGA",
        "GACT", "GCAT", "TACG", "TCAG"))
    t4 <- (!is.null(alias)) && (ncol(mcols(x[[1]])) != length(alias))
    
    if (any(c(t1, t2, t3, t4))) {
        if (t1)
            return(r1)
        if (t2 || t3)
            return(c(r2, r3)[c(t2, t3)])
        if (t4)
            return(message("*** The number of codon sequences must be equal",
                        " to the number of provided sequence's 'alias'."))
    }
    else 
        NULL
}

setValidity2("ListCodonMatrix", valid.ListCodonMatrix)

## ========================== MatrixSeq class =============================

setClassUnion("matrix_OR_vector", c("matrix", "vector", "numeric"))

#' @aliases MatrixSeq
#' @rdname MatrixSeq
#' @title Definition of MatrixSeq-class
#' @description This is a very simple flexible class to store DNA and 
#' aminoacid aligned sequences together with their physicochemical properties.
#' That is, a place where each aminoacid or codon from the sequence is
#' represented by numerical value from a physicochemical index.
#' @details
#' 
#' \describe{
#'  \item{\strong{seqs}: }{A string character vector of DNA or aminoacid 
#'  sequences.}
#'  \item{\strong{matrix}: }{A numerical matrix or a numerical vector
#'  (in the constructor) carrying the specified aminoacid physicochemical
#'  indices for aminoacid in the DNA or aminoacid sequence(s).}
#'  \item{\strong{names}: }{Alias/names/IDs DNA or aminoacid sequences.}
#'  \item{\strong{aaindex}: }{Aminoacid index database where the 
#'  physicochemical index can be found.}
#'  \item{\strong{phychem}: }{Description of the physicochemical index applied 
#'  to represent the DNA or aminoacid sequences.}
#'  \item{\strong{accession}: }{Accession number or ID of the applied
#'  physicochemical index in the database.}
#' }
#' 
#' @keywords internal
#' @author Robersy Sanchez <https://genomaths.com>
#' @export
#' @return Given the slot values, it defines a MatrixSeq-class.
setClass("MatrixSeq",
    slots = c(
        seqs = "character",
        matrix = "matrix",
        names = "character",
        aaindex = "character",
        phychem = "character",
        accession = "character"
    )
)

#' MatrixSeq Constructor
#' @rdname MatrixSeq
#' @param seqs,matrix,names,aaindex,phychem,accession See detail section 
#' @export
#' @returns A MatrixSeq-class object
#' @examples 
#' aln <- c(S1 = "ATGCGGATTAGA", S2 = "ATGACGATCACA", S3 = "ATGAGATCACAG")
#' cd <- DNAMultipleAlignment(aln)
#' r1 <- peptide_phychem_index(unmasked(cd), acc = "EISD840101")
#' r1
#' 
#' ## Extract the second aminoacid sequence
#' r1[2]
#' 
#' ## Using the sequence given name 
#' r1$S1
#' 
#' ## Extract the second aminoacid value from the first sequence
#' r1[1,2]
#' 
#' ## Change the name the second sequence
#' names(r1) <- c('S1', 'Seq1', 'S1')
#' r1
#' 
#' ## Extract the amino acid sequences
#' slot(r1, 'seqs')
#' 
MatrixSeq <- function(seqs, matrix, names, aaindex, phychem, accession) {
    vtr <- is.vector(matrix)
    if (vtr)
        matrix <- as.matrix(matrix)
    
    if (!vtr && nrow(matrix) != length(seqs))
        stop("*** The number sequences must be equal to the matrix ",
            "row-number.")
    
    new("MatrixSeq",
        seqs = seqs,
        matrix = matrix,
        names = names,
        aaindex = aaindex,
        phychem = phychem,
        accession = accession)
}


## ======================= GRangesMatrixSeq class =======================

#' @aliases GRangesMatrixSeq
#' @rdname GRangesMatrixSeq
#' @title Definition of GRangesMatrixSeq-class
#' @description This is a very simple flexible class to store DNA and 
#' aminoacid aligned sequences together with their physicochemical properties.
#' That is, a place where each aminoacid or codon from the sequence is
#' represented by numerical value from a physicochemical index.
#' @keywords internal
#' @export
#' @return Given the slot values, it defines a MatrixList-class.
setClass("GRangesMatrixSeq",
    slots = c(
        seqnames = "Rle",
        ranges = "IRanges_OR_IPos",
        strand = "Rle",
        elementMetadata = "DataFrame",
        seqinfo = "Seqinfo",
        seqs = "character",
        names = "character",
        aaindex = "character",
        phychem = "character",
        accession = "character"
    ),
    contains = "GRanges"
)

#' @importFrom methods new setAs
setAs("MatrixSeq", "GRangesMatrixSeq",
    function(from, to) {

        seqs <- from@seqs
        names <- from@names
        aaindex <- from@aaindex
        phychem <- from@phychem
        accession <- from@accession
        
        pos <- colnames(from@matrix)
        pos <- as.numeric(gsub("A", "", pos))
        seqinfo <- Seqinfo(seqnames = "1", seqlengths = range(pos)[2], 
                        isCircular = NA, genome = NA)
        
        from <- data.frame(
                    seqnames = rep("1", length(pos)), 
                    start = pos, 
                    end = pos, 
                    strand = rep("+", length(pos)),
                    t(from@matrix))
        
        from <- makeGRangesFromDataFrame(from, keep.extra.columns = TRUE)
        
        
        new("GRangesMatrixSeq",
            seqnames = seqnames(from),
            ranges = ranges(from),
            strand = strand(from),
            elementMetadata = from@elementMetadata,
            seqinfo = seqinfo,
            seqs = seqs,
            names = names,
            aaindex = aaindex,
            phychem = phychem,
            accession = accession)
    }
)


#' GRangesMatrixSeq-class constructor
#' @rdname GRangesMatrixSeq
#' @aliases GRangesMatrixSeq
#' @description
#' Constructor for 'GRangesMatrixSeq-class' object.
#' @details
#' This is a convenient function to transform a [MatrixSeq]-class object 
#' returned by function [aa_phychem_index] into a 'GRangesMatrixSeq-class' 
#' object. Since a 'GRangesMatrixSeq-class' inherits from 
#' \code{\link[GenomicRanges]{GRanges-class}}, this transformation permits the 
#' application of several methods from GenomicRanges package in the downstream
#' analysis.
#' @param object If provided, it must be a GRangesMatrixSeq-class object and
#' in this case
#' @param seqnames,start,end,ranges,strand,elementMetadata,seqinfo The same as
#' in \code{\link[GenomicRanges]{GRanges}}
#' @param seqs,names,aaindex,phychem,accession The same as in 
#' [MatrixSeq].
#' @export
#' @examples
#' aln <- c(S1 = "ATGCGGATTAGA", S2 = "ATGACGATCACA", 
#'         S3 = "ATGAGATCACAG")
#' cd <- DNAMultipleAlignment(aln)
#' r1 <- peptide_phychem_index(unmasked(cd), acc = "EISD840101")
#' 
#' r2 <- GRangesMatrixSeq(r1)
#' r2
#' 
#' slot(r2, "phychem")
#' 
GRangesMatrixSeq <- function(
    object = NULL,
    seqnames = Rle(factor()), 
    start = integer(0), 
    end = integer(0),
    ranges = IRanges(),
    strands = Rle(strand()), 
    elementMetadata = DataFrame(),
    seqinfo = NULL,
    seqs = character(),
    names = character(),
    aaindex = character(),
    phychem = character(),
    accession = character()) {
    
    if (!is.null(object)) {
        return(as(object, "GRangesMatrixSeq"))
    }
    else {
        
        if (is.null(seqinfo))
            seqinfo <- seqinfo(ranges)

        new("GRangesMatrixSeq",
            seqnames = seqnames,
            ranges = IRanges(start = start, end = end),
            strand = strands,
            elementMetadata = elementMetadata,
            seqinfo = seqinfo,
            names = names,
            aaindex = aaindex,
            phychem = phychem,
            accession = accession)
    }
}

#' @rdname GRangesMatrixSeq
#' @keywords internal
#' @import Biostrings
#' @export
#' @return Only used to specify signature in the S4 setMethod.
setClassUnion(
    "DNAStringSet_OR_DNAMultipleAlignment",
    c("DNAStringSet", "DNAMultipleAlignment")
)

## ======================== MatrixList class ==========================

#' @aliases MatrixList
#' @rdname MatrixList
#' @title Definition of MatrixList-class
#' @description  A class denoting a list of matrices.
#' @keywords internal
#' @export
#' @return Given the slot values, it defines a MatrixList-class.
setClass("MatrixList",
    slots = c(
        matrices = "list",
        names = "character"
    )
)

## ======================== Validity MatrixList ======================= #
#' @rdname valid.MatrixList
#' @title Valid MatrixList
#' @param x A 'MatrixList object'
#' @keywords internal
#' @return If valid return NULL
valid.MatrixList <- function(x) {
    if (!all(slapply(x, function(y) inherits(y, "matrix")))) {
        return(
            "*** Not all the elements of the MatrixList object",
            " are from 'matrix' class."
        )
    }
    NULL
}


#' @rdname valid.MatrixList
#' @keywords internal
#' @import Biostrings
#' @export
#' @return Only used to specify signature in the S4 setMethod.
setClassUnion(
    "DNAStringSet_OR_NULL",
    c("DNAStringSet", "DNAMultipleAlignment", "NULL", "missing")
)


## ========================== Automorphism =============================

#' @rdname Automorphism
#' @aliases Automorphism
#' @title A class definition to store codon automorphisms in a given Abelian
#' group representation.
#' @description Two classes are involved in to storing codon automorphisms:
#' \emph{\strong{Automorphism-class}} and
#' \emph{\strong{AutomorphismList-class}}.
#' @details An \emph{\strong{Automorphism-class}} object has six 
#' columns: "seq1", "seq2","coord1", "coord2", "autm", and "cube". See the
#' examples for function \code{\link{automorphisms}}. Observe that as the
#' \emph{\strong{Automorphism-class}} inherits from
#' \code{\link[GenomicRanges]{GRanges-class}} the transformation starting 
#' from a \code{\link[GenomicRanges]{GRanges-class}} object into an
#' \emph{\strong{Automorphism-class}} is straightforward. 
#' 
#' However, the transformation starting from a \code{\link[base]{data.frame}}
#' or a \code{\link[S4Vectors]{DataFrame-class}} object \eqn{"x"} requires for
#' the creation of an additional \code{\link[GenomicRanges]{GRanges-class}}
#' object, which by default will have the argument seqnames = "1", strand =
#' "+", start/end = seq(row(x)), length = nrow(x). These details must be keep
#' in mind to prevent fundamental errors in the downstream analyses.
#' 
#' @section Automorphism-class methods:
#' ## as(from, "Automorphism"):
#' Permits the transformation of a \code{\link[base]{data.frame}} or a
#' \code{\link[S4Vectors]{DataFrame-class}} object into
#' \emph{\strong{Automorphism-class}} object if the proper columns are 
#' provided. 
#' 
#' Methods from \code{\link[GenomicRanges]{GRanges-class}} can also be 
#' applied.
#' @keywords internal
#' @import S4Vectors
#' @import GenomicRanges
#' @seealso \code{\link{AutomorphismByCoef-class}} and
#' \code{\link{AutomorphismList-class}}
#' @export
#' @return Given the slot values, it defines an Automorphism-class object.
setClass("Automorphism",
    slots = c(
        seqnames = "Rle",
        ranges = "IRanges_OR_IPos",
        strand = "Rle",
        elementMetadata = "DataFrame",
        seqinfo = "Seqinfo",
        colnames = "character",
        autm_info = "list"
    ),
    contains = "GRanges"
)

#' @rdname Automorphism
#' @import S4Vectors
#' @keywords internal
#' @export
setClassUnion(
    "DataFrame_OR_data.frame",
    c("DataFrame", "data.frame")
)


#' @importFrom GenomeInfoDb Seqinfo seqnames
#' @import S4Vectors
#' @import GenomicRanges
#' @import Biostrings
#' @importFrom BiocGenerics strand
#' @importFrom IRanges IRanges ranges
#' @importFrom methods new setAs
setAs(
    "DataFrame_OR_data.frame", "Automorphism",
    function(from) {
        nr <- nrow(from)
        pos <- seq(1, nr, 1)
        gr <- GRanges(
            seqnames = 1,
            ranges = IRanges(start = pos, end = pos),
            strand = "+"
        )
        mcols(gr) <- from

        new("Automorphism",
            seqnames = seqnames(gr),
            ranges = ranges(gr),
            strand = strand(gr),
            elementMetadata = from,
            seqinfo = seqinfo(gr),
            colnames = colnames(from)
        )
    }
)


# ======================== Validity Automorphism ======================= #
#' @rdname valid.Automorphism
#' @title Valid Automorphism mcols
#' @param x A 'Automorphism object'
#' @keywords internal
#' @return An Error if the metacolumn does not have a valid format
valid.Automorphism.mcols <- function(x) {
    alf <- c("A", "C", "G", "T", "-")
    if (length(x) > 0) {
        m1 <- m2 <- FALSE
        if (inherits(x, "GRanges")) {
            coln <- colnames(x@elementMetadata)
        } else {
            coln <- x@colnames
        }
        if (any(!is.element(
            c(
                "seq1", "seq2", "coord1",
                "coord2", "autm", "cube"
            ),
            coln
        ))) {
            m1 <- TRUE
        }
        if (unique(nchar(x$seq1)) != 3 || unique(nchar(x$seq2)) != 3) {
            m2 <- TRUE
        }
        if (m2) {
            if (all(is.element(x$seq1, alf)) &&
                all(is.element(x$seq1, alf))) {
                m2 <- FALSE
            }
        }

        if (m1 || m2) {
            return("*** This is not a valid Automorphism-class object.")
        }
    }
    NULL
}

#' @rdname valid.Automorphism
#' @title Valid Automorphism
#' @param x A 'Automorphism object'
#' @keywords internal
#' @return An Error if the Automorphism-class object is not valid.
valid.Automorphism <- function(x) {
    c(valid.GRanges(x), valid.Automorphism.mcols(x))
}

setValidity2("Automorphism", valid.Automorphism)


## ======================= AutomorphismList-class =========================

#' @rdname AutomorphismList
#' @title A class definition to store list of Automorphism class objects.
#' @description A class definition to store list of Automorphism class objects
#' derived from the pairwise automorphism estimation from pairwise
#' alignments. Objects from this class are created by function 
#' \code{\link{automorphisms}} and \code{\link{as.AutomorphismList}}.
#' @importFrom methods validObject setClass
#' @keywords internal
#' @return An object from AutomorphismList-class 
#' @export
#' @aliases AutomorphismList
#' @section AutomorphismList-class methods:
#' ## as.AutomorphismList(x):
#' \emph{\strong{as.AutomorphismList}} function transform a list of
#' \code{\link[GenomicRanges]{GRanges-class}}, a
#' \code{\link[GenomicRanges]{GRangesList-class}}, a list of
#' \code{\link[base]{data.frame}} or a
#' \code{\link[S4Vectors]{DataFrame-class}}
#' objects into a \emph{\strong{AutomorphismList-class}} object.
#' 
#' ## unlist(x)
#' It transforms a AutomorphismList-class object into an Automorphism-class
#' object. 
#'
#' ## as.list(x)
#' It transforms a list of Automorphism-class objects into an 
#' AutomorphismList-class object.
#' 
#' ## as(x, "GRangesList")
#' It transforms a \code{\link[GenomicRanges]{GRangesList}} of 
#' \code{\link{Automorphism-class}} objects into an 
#' 'AutomorphismList-class' object.
#
#' ## names(x)
#' To get the element's names from an 'AutomorphismList-class' object.
#'
#' ## names(x) <- value
#' To assign names to the element from an 'AutomorphismList-class' 
#' object.
#' @examples
#' ## Load datasets
#' data("autm", "brca1_autm")
#' 
#' ## Transforming a list of Automorphisms into an AutomorphismList object
#' lista <- list(human = brca1_autm[[1]], gorilla = brca1_autm[[2]])
#' as.AutomorphismList(lista)
#' 
#' ## Alternatively we can set
#' aut <- as.list(brca1_autm[seq(2)])
#' class(aut)
#' 
#' ## And reverse it
#' aut <- as.AutomorphismList(aut)
#' aut
#' 
#' ## Let's get the element names from object 'aut'
#' names(aut)
#' 
#' ## Let's assign new names
#' names(aut) <- c("human_1", "gorilla_1")
#' names(aut)
#' 
#' ## Transforming a GRangesList of Automorphisms into an AutomorphismList
#' ## object
#' lista <- as(lista, "GRangesList")
#' as.AutomorphismList(lista)
#'
#' ## Transform a AutomorphismList-class object into an Automorphism-class
#' ## object 
#' unlist(brca1_autm[seq(2)])
#' @seealso \code{\link{Automorphism-class}} and 
#' \code{\link{AutomorphismByCoefList-class}}.
setClass("AutomorphismList",
    slots = c(
        DataList = "list",
        SeqRanges = "GenomicRanges_OR_missing"
    )
)

## ====================== Validity AutomorphismList ================== #

#' @rdname valid.AutomorphismList
#' @title Valid AutomorphismList mcols
#' @param x A 'AutomorphismList object'
#' @return An error if 'x' is not a valid AutomorphismList class object.
#' @import S4Vectors
#' @keywords internal

valid.AutomorphismList <- function(x) {
    m1 <- FALSE
    if (!(inherits(x@DataList[[1]], "Automorphism") ||
        inherits(x@DataList[[1]], "DataFrame"))) {
        m1 <- TRUE
    }

    if (inherits(x@DataList[[1]], "Automorphism")) {
        if (any(!slapply(
            x@DataList,
            function(y) {
                return(inherits(y, "Automorphism") && validObject(y))
            }
        ))) {
            m1 <- TRUE
        }
    }

    if (inherits(x@DataList[[1]], "DataFrame")) {
        if (any(!slapply(
            x@DataList,
            function(y) {
                return(inherits(y, "DataFrame") && validObject(y))
            }
        ))) {
            m1 <- TRUE
        }
        if (!inherits(x@SeqRanges, "GRanges")) {
            m1 <- TRUE
        }
    }
    if (m1) {
        return("*** This is not a valid AutomorphismList class object.")
    }
    NULL
}

setValidity2("AutomorphismList", valid.AutomorphismList)

#' @rdname AutomorphismList
#' @param x An \code{\link{AutomorphismList}} object.
#' @export
#' @examples 
#' ## Load a DNA sequence alignment
#' data("brca1_autm", package = "GenomAutomorphism")
#' names(brca1_autm)
setMethod("names",
        signature = "AutomorphismList",
    function(x) names(x@DataList)
)

## ================= AutomorphismList-methods ========================

#' @rdname AutomorphismList
#' @aliases 'names<-'
#' @param x An \code{\link{AutomorphismList-class}} object.
#' @param value A character vector naming the elements of the 
#' \code{\link{AutomorphismList-class}} object 'x'.
#' @export
#' @examples 
#' ## Load a DNA sequence alignment
#' data("brca1_autm", package = "GenomAutomorphism")
#' x1 <- brca1_autm[seq(2)]
#' names(x1)
#' 
#' ## Let's assign a new names
#' names(x1) <- c("human_1.human_2.0", "human_1.gorilla_0")
#' names(x1) 
setReplaceMethod(
    "names", "AutomorphismList",
    function(x, value) {
        names(x@DataList) <- value
        return(x)
    }
)

#' @rdname AutomorphismList
#' @export
#' @examples 
#' ## Load a DNA sequence alignment
#' data("brca1_autm", package = "GenomAutomorphism")
#' 
#' ## The list of the first three elements
#' autm_list <- as.list(brca1_autm[seq(3)])
#' autm_list
setMethod("as.list",
    signature = "AutomorphismList",
    function(x) {
        x <- getAutomorphisms(x)
        return(x@DataList)
    }
)


#' @importFrom methods setAs coerce
setAs(from = "AutomorphismList", to = "list", function(from) {
    from <- getAutomorphisms(from)
    return(from@DataList)
})

#' @importFrom methods setAs
setAs("AutomorphismList", "GRangesList", function(from) {
    from <- getAutomorphisms(from)
    from <- as.list(from)
    return(as(from, "GRangesList"))
})


#' @import GenomicRanges
setMethod("unlist",
    signature = "AutomorphismList",
    function(x) {
        x <- as(x, "GRangesList")
        return(unlist(x))
    }
)

## ========================== AutomorphismByCoef ===========================

#' @aliases AutomorphismByCoef
#' @rdname AutomorphismByCoef
#' @title A class definition to store conserved gene/genomic regions found
#' in a MSA. 
#' @description Objects from this class are generated by function 
#' \code{\link{automorphism_bycoef}}.
#' @seealso \code{\link{automorphism_bycoef}}
#' @keywords internal
#' @import GenomicRanges
#' @section AutomorphismByCoefList-class methods:
#' ## unlist(x):
#' It transforms a AutomorphismByCoefList-class object into an 
#' AutomorphismByCoef-class object. 
#' 
#' ## as(x, "AutomorphismByCoefList")
#' It transforms a 'list' of AutomorphismByCoef-class object into an 
#' AutomorphismByCoefList-class object.
#' @export
#' @examples 
#' ## Let's transform a AutomorphismByCoefList-class object into an 
#' ## AutomorphismByCoef-class object
#' data("autby_coef")
#' unlist(autby_coef[1:2])
#' 
#' ## Herein a 'list' object of AutomorphismByCoef-class objects
#' lista <- list(human = autby_coef[[1]], gorilla = autby_coef[[2]])
#' 
#' ## Let's transform the the last list 'lista' into an
#' ## AutomorphismByCoefList-class object
#' aut <- as(lista, "AutomorphismByCoefList")
#' aut
#' 
#' ## Let's get the element names from object 'aut'
#' names(aut)
#' 
#' ## Let's assign new names
#' names(aut) <- c("human_1", "gorilla_1")
#' names(aut)
#' @seealso \code{\link{AutomorphismByCoefList-class}} and 
#' \code{\link{Automorphism-class}}
#' @return AutomorphismByCoef-class definition.
setClass("AutomorphismByCoef",
    slots = c(
        seqnames = "Rle",
        ranges = "IRanges_OR_IPos",
        strand = "Rle",
        elementMetadata = "DataFrame",
        seqinfo = "Seqinfo",
        colnames = "character",
        autm_info = "list"),
    contains = "GRanges"
)

# ======================== Validity AutomorphismByCoef ================== #
#' @rdname valid.AutomorphismByCoef
#' @aliases valid.AutomorphismByCoef
#' @title Valid AutomorphismByCoef mcols
#' @param x A 'AutomorphismByCoef object'
#' @import S4Vectors
#' @return An error if 'x' is not a valid AutomorphismByCoef.
#' @keywords internal
valid.AutomorphismByCoef <- function(x) {
    coln <- colnames(mcols(x))
    if (!inherits(x, "GRanges") ||
        any(!is.element(c("autm", "cube"), coln))) {
        return("*** This is not a valid AutomorphismByCoef
                class object.")
    } else {
        NULL
    }
}

setValidity2("AutomorphismByCoef", valid.AutomorphismByCoef)


## ========================= AutomorphismByCoefList ======================

#' @aliases AutomorphismByCoefList
#' @rdname AutomorphismByCoefList
#' @title A class definition for a list of AutomorphismByCoef class objects.
#' @keywords internal
#' @import S4Vectors
#' @importFrom methods as
#' @details \strong{AutomorphismByCoefList-class} has the following methods:
#' ## as('from', "AutomorphismByCoefList")
#' Where 'from' is a list of \strong{AutomorphismByCoef-class}.
#'
#' ## unlist(x)
#' Where 'x' is a an \strong{AutomorphismByCoefList-class} object.
#' @export
#' @seealso \code{\link{AutomorphismByCoef-class}} and 
#' \code{\link{AutomorphismList-class}}
#' @return AutomorphismByCoefList-class definition.
setClass(
    "AutomorphismByCoefList",
    slots = c(
        elementMetadata = "DataFrame",
        elementType = "character",
        metadata = "list",
        listData = "list"
    ),
    contains = "SimpleGRangesList"
)

# ===================== Validity AutomorphismByCoefList ================== #
#' @rdname valid.AutomorphismByCoefList
#' @aliases valid.AutomorphismByCoefList
#' @title Valid AutomorphismByCoefList mcols
#' @param x A 'AutomorphismByCoefList object'
#' @import S4Vectors
#' @keywords internal
#' @return An error if 'x' is not a valid AutomorphismByCoefList.
valid.AutomorphismByCoefList <- function(x) {
    if (any(!slapply(x, validObject)) || any(slapply(x, function(y) {
        coln <- colnames(mcols(y))
        !is.element(c("autm", "cube"), coln)
    }))) {
        return("*** This is not a valid AutomorphismByCoefList
                class object.")
    } else {
        NULL
    }
}

setValidity2(
    "AutomorphismByCoefList",
    valid.AutomorphismByCoefList
)


as_list_of_AutomorphismByCoef <- function(from) {
    lapply(from, as, Class = "AutomorphismByCoef")
}

#' @importFrom methods setAs
setAs("list", "AutomorphismByCoefList", function(from) {
    from <- as_list_of_AutomorphismByCoef(from)
    from <- new_SimpleList_from_list(
        Class = "SimpleGRangesList",
        x = from
    )
    new("AutomorphismByCoefList", from)
})

setMethod("unlist",
    signature = "AutomorphismByCoefList",
    function(x) {
        x <- as(x, "GRangesList")
        return(unlist(x))
    }
)

## ========================== ConservedRegion ===========================

#' @aliases ConservedRegion
#' @rdname ConservedRegion
#' @title A class definition to store conserved gene/genomic regions found
#' in a MSA.
#' @keywords internal
#' @export
#' @return Definition of the \strong{ConservedRegion-class}.
setClass("ConservedRegion",
    contains = "GRanges"
)

# ======================== Validity ConservedRegion ================== #
#' @aliases valid.ConservedRegion
#' @rdname ConservedRegion
#' @title Valid ConservedRegion mcols
#' @param x A 'ConservedRegion object'
#' @import S4Vectors
#' @keywords internal

valid.ConservedRegion <- function(x) {
    coln <- colnames(mcols(x))
    if (!inherits(x, "GRanges") ||
        any(!is.element(coln, c("autm", "cube")))) {
        return("*** This is not a valid ConservedRegion
                class object.")
    } else {
        NULL
    }
}


setValidity2("ConservedRegion", valid.ConservedRegion)

## ========================= ConservedRegionList ======================

#' @aliases ConservedRegionList
#' @rdname ConservedRegion
#' @title A class definition for a list of ConservedRegion class objects.
#' @details \strong{ConservedRegionList-class} has the following method:
#' ## as('from', "ConservedRegionList")
#' Where 'from' is a list of \strong{ConservedRegion-class}.
#'
#' @keywords internal
#' @export
setClass(
    "ConservedRegionList",
    slots = c(
        elementMetadata = "DataFrame",
        elementType = "character",
        metadata = "list",
        listData = "list"
    ),
    contains = "SimpleGRangesList"
)

as_list_of_ConservedRegion <- function(from) {
    lapply(from, as, Class = "ConservedRegion")
}

setAs("list", "ConservedRegionList", function(from) {
    from <- as_list_of_ConservedRegion(from)
    from <- new_SimpleList_from_list(
        Class = "SimpleGRangesList",
        x = from
    )
    new("ConservedRegionList", from)
})

# ===================== Validity ConservedRegionList ================== #
#' @rdname ConservedRegion
#' @aliases valid.ConservedRegion
#' @title Valid ConservedRegionList mcols
#' @param x A 'ConservedRegionList object'
#' @import S4Vectors
#' @keywords internal

valid.ConservedRegionList <- function(x) {
    coln <- colnames(mcols(x))
    if (any(!slapply(x, validObject)) || any(slapply(x, function(y) {
        coln != "autm"
    }))) {
        return("*** This is not a valid ConservedRegionList
                class object.")
    } else {
        NULL
    }
}

setValidity2("ConservedRegionList", valid.ConservedRegionList)

## ========================== Class union ============================ 

setClassUnion("CodonGroup_OR_Automorphisms", 
            c("CodonGroup", "Automorphism", "AutomorphismByCoef"))

