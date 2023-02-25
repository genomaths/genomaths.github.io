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

#' @rdname codon_coord
#' @title Codon coordinates on a given a given Abelian group representation.
#' @description Given a string denoting a codon or base from the DNA (or RNA)
#' alphabet and a genetic-code Abelian group as given in reference (1).
#' @param codon An object from \code{\link{BaseGroup-class}} (generated with
#' function \code{\link{base_coord}}), \code{\link[Biostrings]{DNAStringSet}}
#' or from \code{\link[Biostrings]{DNAMultipleAlignment}} class carrying the
#' DNA pairwise alignment of two sequences.
#' @param filepath A character vector containing the path to a file in
#' \emph{\strong{fasta}} format to be read. This argument must be given if
#' \emph{codon & base} arguments are not provided.
#' @param cube A character string denoting one of the 24 Genetic-code cubes,
#' as given in references (2-3).
#' @param group A character string denoting the group representation for the
#' given base or codon as shown in reference (2-3).
#' @param start,end,chr,strand Optional parameters required to build a
#' \code{\link[GenomicRanges]{GRanges-class}}. If not provided the default
#' values given for the function definition will be used.
#' @param ... Not in use.
#' @details Symbols "-" and "N" usually found in DNA sequence alignments to
#' denote gaps and missing/unknown bases are represented by the number: '-1'
#' on Z4 and '0' on Z5. In Z64 the symbol 'NA' will be returned for codons
#' including symbols "-" and "N".
#'
#' This function returns a \code{\link[GenomicRanges]{GRanges-class}} object
#' carrying the codon sequence(s) and their respective coordinates in the
#' requested Abelian group or simply, when \emph{group =  'Z5^3'}
#' 3D-coordinates, which are derive from Z5 as indicated in reference (3).
#' Notice that the coordinates can be 3D or just one-dimension ("Z64" or
#' "Z125"). Hence, the pairwise alignment provided in argument
#' \emph{\strong{codon}} must correspond to codon sequences.
#'
#' @seealso [Symmetric Group of the Genetic-Code Cubes.](
#' https://github.com/genomaths/GenomeAlgebra_SymmetricGroup)
#' @import S4Vectors
#' @import Biostrings
#' @importFrom methods new
#' @export
#' @return A \code{\link{CodonGroup-class}} object.
#' @seealso \code{\link{base_coord}} and \code{\link{base2int}}.
#' @author Robersy Sanchez <https://genomaths.com>
#' @references
#' \enumerate{
#'  \item Robersy Sanchez, Jesus Barreto (2021) Genomic Abelian Finite
#'   Groups.
#'  [doi: 10.1101/2021.06.01.446543](https://doi.org/10.1101/2021.06.01.446543)
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
#' ## DNA base representation in the Abelian group Z5
#' bs_cor <- codon_coord(
#'     codon = aln,
#'     cube = "ACGT",
#'     group = "Z5"
#' )
#' bs_cor ## 3-D coordinates
#'
#'
#' ## DNA base representation in the Abelian group Z64
#' bs_cor <- codon_coord(
#'     codon = aln,
#'     cube = "ACGT",
#'     group = "Z64"
#' )
#' bs_cor
#'
#' ## Giving a matrix of codons
#' codon_coord(base2codon(x = aln))
#'
#' @aliases codon_coord
setGeneric(
    "codon_coord",
    function(codon = NULL,
    ...) {
        standardGeneric("codon_coord")
    }
)

#' @aliases codon_coord
#' @rdname codon_coord
#' @import GenomicRanges
#' @import BiocGenerics
#' @importFrom GenomeInfoDb seqnames
#' @import IRanges
#' @import S4Vectors
setMethod(
    "codon_coord", signature(codon = "BaseGroup"),
    function(codon,
    group = NULL) {
        if (is.null(group)) {
            group <- codon@group
        }
        cube <- codon@cube
        chr <- unique(as.character(seqnames(codon)))
        strands <- unique(as.character(strand(codon)))

        codon <- mcols(codon)
        idx_seq <- grep("seq", colnames(codon))
        idx_coord <- grep("coord", colnames(codon))

        sq <- data.frame(codon[, idx_seq])
        f <- factor(as.vector(slapply(seq_len(nrow(sq) / 3), rep, times = 3)))
        sq <- split(sq, f)

        crd <- data.frame(codon[, idx_coord])
        crd <- split(crd, f)

        idx <- seq_along(sq)
        if (is.element(group, c("Z4", "Z5", "Z4^3", "Z5^3"))) {
            codon <- lapply(idx, function(k) {
                c(
                    apply(sq[[k]], 2, paste, collapse = ""),
                    apply(crd[[k]], 2, paste, collapse = ",")
                )
            })
        } else {
            fun <- switch(group,
                Z64 = CodonCoordZ4toZ64,
                Z125 = CodonCoordZ5toZ125,
            )

            codon <- lapply(idx, function(k) {
                c(
                    apply(sq[[k]], 2, paste, collapse = ""),
                    apply(crd[[k]], 2, fun)
                )
            })
        }

        rm(sq, crd)
        gc()
        codon <- do.call(rbind, codon)

        pos <- seq(1, nrow(codon), 1L)
        codon <- data.frame(
            seqnames = chr,
            start = pos,
            end = pos,
            strand = strands,
            codon
        )

        codon <- makeGRangesFromDataFrame(codon, keep.extra.columns = TRUE)
        codon <- new2(
            "CodonGroup",
            seqnames = seqnames(codon),
            ranges = ranges(codon),
            strand = strand(codon),
            elementMetadata = codon@elementMetadata,
            seqinfo = codon@seqinfo,
            colnames = colnames(codon@elementMetadata),
            group = group,
            cube = cube, check = FALSE
        )
        return(codon)
    }
)


#' @aliases codon_coord
#' @rdname codon_coord
#' @import Biostrings
setMethod(
    "codon_coord", signature(codon = "DNAStringSet_OR_NULL"),
    function(codon = NULL,
    filepath = NULL,
    cube = c(
        "ACGT", "AGCT", "TCGA", "TGCA", "CATG",
        "GTAC", "CTAG", "GATC", "ACTG", "ATCG",
        "GTCA", "GCTA", "CAGT", "TAGC", "TGAC",
        "CGAT", "AGTC", "ATGC", "CGTA", "CTGA",
        "GACT", "GCAT", "TACG", "TCAG"
    ),
    group = c("Z4", "Z5", "Z64", "Z125", "Z4^3", "Z5^3"),
    start = NA,
    end = NA,
    chr = 1L,
    strand = "+") {
        if (!is.null(filepath) && is.character(filepath)) {
            codon <- readDNAMultipleAlignment(filepath = filepath)
        }

        if (any(nchar(codon) %% 3 != 0)) {
            stop(
                "*** 'codon' argument is not a base-triplet sequence.",
                " A base-triplet sequence is multiple of 3."
            )
        }

        cube <- match.arg(cube)
        group <- match.arg(group)

        base_grp <- group
        if (is.element(group, "Z64")) {
            base_grp <- "Z4"
        }
        if (is.element(group, "Z125")) {
            base_grp <- "Z5"
        }
        if (is.element(group, "Z4^3")) {
            base_grp <- "Z4"
        }
        if (is.element(group, "Z5^3")) {
            base_grp <- "Z5"
        }

        codon <- base_coord(
            base = codon,
            filepath = NULL,
            cube = cube,
            group = base_grp,
            start = start,
            end = end,
            chr = chr,
            strand = strand
        )
        codon <- codon_coord(codon, group = group)
        return(codon)
    }
)

setClassUnion("matrix_OR_data_frame", c("matrix", "data.frame"))

#' @aliases codon_coord
#' @rdname codon_coord
setMethod(
    "codon_coord", signature(codon = "matrix_OR_data_frame"),
    function(codon,
        cube = c(
            "ACGT", "AGCT", "TCGA", "TGCA", "CATG",
            "GTAC", "CTAG", "GATC", "ACTG", "ATCG",
            "GTCA", "GCTA", "CAGT", "TAGC", "TGAC",
            "CGAT", "AGTC", "ATGC", "CGTA", "CTGA",
            "GACT", "GCAT", "TACG", "TCAG"
        ),
        group = c("Z64", "Z125", "Z4^3", "Z5^3")) {
        
        cube <- match.arg(cube)
        group <- match.arg(group)

        crd <- data.frame(codon)
        base_grp <- group
        if (is.element(group, "Z64")) {
            base_grp <- "Z4"
        }
        if (is.element(group, "Z125")) {
            base_grp <- "Z5"
        }
        if (is.element(group, "Z4^3")) {
            base_grp <- "Z4"
        }
        if (is.element(group, "Z5^3")) {
            base_grp <- "Z5"
        }

        crd <- lapply(seq_len(nrow(codon)), function(k) {
            c1 <- base2int(
                base = crd[k, 1],
                cube = cube,
                group = base_grp
            )

            c2 <- base2int(
                base = crd[k, 2],
                cube = cube,
                group = base_grp
            )

            if (is.element(group, c("Z64", "Z125"))) {
                c1 <- switch(group,
                    "Z64" = CodonCoordZ4toZ64(c1),
                    "Z125" = CodonCoordZ5toZ125(c1)
                )
                c2 <- switch(group,
                    "Z64" = CodonCoordZ4toZ64(c2),
                    "Z125" = CodonCoordZ5toZ125(c2)
                )
            } else {
                c1 <- paste(c1, collapse = ",")
                c2 <- paste(c2, collapse = ",")
            }
            return(c(c1, c2))
        })
        crd <- do.call(rbind, crd)

        crd <- data.frame(codon, crd)
        colnames(crd) <- c(paste0("seq", seq(dim(codon))),
                            paste0("coord", seq(dim(codon))))
        return(crd)
    }
)


## ------------------------- Auxiliary functions --------------------------

CodonCoordZ4toZ64 <- function(x) {
    if (any(is.na(x))) {
        res <- NA
    } else {
        res <- 4 * x[1] + 16 * x[2] + x[3]
    }
    return(res)
}

CodonCoordZ5toZ125 <- function(x) {
    if (any(is.na(x))) {
        res <- NA
    } else {
        res <- 5 * x[1] + 25 * x[2] + x[3]
    }
    return(res)
}
