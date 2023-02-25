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

#' @rdname autZ64
#' @aliases autZ64
#' @title Compute the Automorphisms of Mutational Events Between two Codon
#' Sequences Represented in Z64.
#' @description Given two codon sequences represented in the Z64 Abelian group,
#' this function computes the automorphisms describing codon mutational events.
#' @details Automorphisms in Z64 are described as functions
#' \eqn{f(x) = k * x} \emph{mod 64}, where \eqn{k} and \eqn{x} are elements 
#' from the set of integers modulo 64.
#' 
#' @param seq An object from a \code{\link[Biostrings]{DNAStringSet}} or
#' \code{\link[Biostrings]{DNAMultipleAlignment}} class carrying the DNA
#' pairwise alignment of two sequences. The pairwise alignment provided in
#' argument \emph{\strong{seq}} or the 'fasta' file \emph{\strong{filepath}}
#' must correspond to codon sequences.
#' @param filepath A character vector containing the path to a file in
#' \emph{\strong{fasta}} format to be read. This argument must be given if
#' \emph{codon & base} arguments are not provided.
#' @param cube,cube_alt A character string denoting pairs of the 24
#' Genetic-code cubes, as given in references (2-3). That is, the base pairs
#' from the given cubes must be complementary each other. Such a cube pair are
#' call dual cubes and, as shown in reference (3), each pair integrates group.
#' @param start,end,chr,strand Optional parameters required to build a
#' \code{\link[GenomicRanges]{GRanges-class}}. If not provided the default
#' values given for the function definition will be used.
#' @param genetic_code The named character vector returned by  
#' \code{\link[Biostrings]{getGeneticCode}} or similar. The translation of
#' codon into aminoacids is a valuable information useful for downstream
#' statistical analysis. The standard genetic code is the default argument
#' value applied in the translation of codons into aminoacids
#' (see \code{\link[Biostrings]{GENETIC_CODE_TABLE}}. 
#' @param num.cores,tasks Parameters for parallel computation using package
#' \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#' use, i.e. at most how many child processes will be run simultaneously (see
#' \code{\link[BiocParallel]{bplapply}} and the number of tasks per job (only
#' for Linux OS).
#' @param verbose If TRUE, prints the progress bar.
#' @return An object \code{\link{Automorphism-class}} with four columns on its
#' metacolumn named: \emph{seq1}, \emph{seq2}, \emph{autm}, and \emph{cube}.
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @import Biostrings
#' @importFrom methods new
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @references
#' \enumerate{
#'  \item Sanchez R, Morgado E, Grau R. Gene algebra from a genetic code
#'  algebraic structure. J Math Biol. 2005 Oct;51(4):431-57.
#'  doi: 10.1007/s00285-005-0332-8. Epub 2005 Jul 13. PMID: 16012800. (
#'  [PDF](https://arxiv.org/pdf/q-bio/0412033.pdf)).
#'  \item Robersy Sanchez, Jesus Barreto (2021) Genomic Abelian Finite
#'   Groups.
#'  [doi:10.1101/2021.06.01.446543](https://doi.org/10.1101/2021.06.01.446543)
#'  \item M. V Jose, E.R. Morgado, R. Sanchez, T. Govezensky, The 24 possible
#'  algebraic representations of the standard genetic code in six or in three
#'  dimensions, Adv. Stud. Biol. 4 (2012) 110-152.[PDF](https://is.gd/na9eap).
#'  \item R. Sanchez. Symmetric Group of the Genetic-Code Cubes. Effect of the
#'  Genetic-Code Architecture on the Evolutionary Process MATCH Commun. Math.
#'  Comput. Chem. 79 (2018) 527-560. [PDF](https://bit.ly/2Z9mjM7)
#' }
#' @examples
#' ## Load a pairwise alignment
#' data(aln, package = "GenomAutomorphism")
#' aln
#'
#' ## Automorphism on Z64
#' autms <- autZ64(seq = aln, verbose = FALSE)
#' autms
#' 
autZ64 <- function(seq = NULL,
    filepath = NULL,
    cube = c("ACGT", "TGCA"),
    cube_alt = c("CATG", "GTAC"),
    start = NA,
    end = NA,
    chr = 1L,
    strand = "+",
    genetic_code = getGeneticCode("1"),
    num.cores = detectCores() - 1,
    tasks = 0L,
    verbose = TRUE) {
    if (is.null(filepath) && is.null(seq)) {
        stop("*** One of the arguments 'seq' or 'filepath' must be given.")
    }
    
    if (!is.null(filepath) && is.character(filepath)) {
        seq <- readDNAMultipleAlignment(filepath = filepath)
    }

    if (!is.null(seq)) {
        if (!inherits(seq, c("DNAStringSet", "DNAMultipleAlignment"))) {
            stop(
                "*** Agument 'seq' must belong to 'DNAStringSet'",
                " DNAMultipleAlignment class."
            )
        }

        if (inherits(seq, "DNAMultipleAlignment")) {
            seq <- unmasked(seq)
        }

        if (any(nchar(seq) %% 3 != 0)) {
            stop(
                "*** The argument of 'seq' must be a pairwise alignment",
                " of codon sequences."
            )
        }
    }

    autm1 <- automorfismos(
        seq = seq,
        filepath = NULL,
        cube = cube,
        start = start,
        end = end,
        chr = chr,
        strand = strand,
        genetic_code = genetic_code,
        num.cores = num.cores,
        tasks = tasks,
        verbose = verbose
    )

    idx <- which(is.na(autm1$autm))

    if (length(idx) > 0) {
        autm2 <- automorfismos(
            seq = seq,
            filepath = NULL,
            cube = cube_alt,
            start = start,
            end = end,
            chr = chr,
            strand = strand,
            genetic_code = genetic_code,
            num.cores = num.cores,
            tasks = tasks,
            verbose = verbose
        )
        autm1[idx, ] <- autm2[idx, ]
    }
    autm1 <- new(
        "Automorphism",
        seqnames = seqnames(autm1),
        ranges = ranges(autm1),
        strand = strand(autm1),
        elementMetadata = autm1@elementMetadata,
        seqinfo = autm1@seqinfo,
        colnames = colnames(autm1@elementMetadata),
        autm_info = list(
            cube = cube,
            cube_alt = cube_alt,
            genetic_code = genetic_code)
    )
    return(autm1)
}



## ===================== Auxiliary functions ===========================

#' @importFrom BiocGenerics strand<-
automorfismos <- function(seq,
    filepath,
    cube,
    start = NA,
    end = NA,
    chr = 1L,
    strand = "+",
    genetic_code = getGeneticCode("1"),
    num.cores,
    tasks,
    verbose) {
    sq <- get_coord(
        x = seq,
        output = "all",
        base_seq = FALSE,
        filepath = filepath,
        cube = cube[1],
        group = "Z64",
        start = start,
        end = end,
        chr = chr,
        strand = strand
    )
    
    sq_alt <- get_coord(
        x = seq,
        output = "all",
        base_seq = FALSE,
        filepath = filepath,
        cube = cube[2],
        group = "Z64",
        start = start,
        end = end,
        chr = chr,
        strand = strand
    )
    
    gr <- sq@SeqRanges
    gr$aa1 <- translation(gr$seq1, genetic.code = genetic_code)
    gr$aa2 <- translation(gr$seq2, genetic.code = genetic_code)
    gr$coord1 <- sq@CoordList$coord1
    gr$coord2 <- sq@CoordList$coord2
    gr$autm <- 1
    gr$cube <- cube[1]
    strands <- as.character(strand(gr))

    idx <- sq@CoordList$coord1 != sq@CoordList$coord2
    idx <- sort(c(which(idx), which(is.na(idx))))

    # ## -------------- Setting parallel computation ----------------- #

    progressbar <- FALSE
    if (verbose) {
        progressbar <- TRUE
    }
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

    # ## -------------------------------------------------------------- #


    if (length(idx) != 0) {
        sq <- bplapply(idx, function(k) {
            c1 <- sq@CoordList$coord1[k]
            c2 <- sq@CoordList$coord2[k]

            if (is.na(c1) && is.na(c2)) {
                s <- c(-1, "Gaps")
            } else {
                s <- try(modeq(c1, c2, 64)[1],
                    silent = TRUE
                )

                if (any(s == -1) || inherits(s, "try-error")) {
                    sq <- sq_alt

                    c1 <- sq@CoordList$coord1[k]
                    c2 <- sq@CoordList$coord2[k]

                    s <- try(modeq(c1, c2, 64)[1],
                        silent = TRUE
                    )

                    if (s != -1 && !inherits(s, "try-error")) {
                        s <- c(s, cube[2])
                    }
                } else {
                    s <- c(s, cube[1])
                }
                if (any(s == -1) || inherits(s, "try-error")) {
                    s <- c(0, "Trnl")
                }
            }
            return(s)
        }, BPPARAM = bpparam)

        sq <- do.call(rbind, sq)
        sq <- data.frame(sq)
        colnames(sq) <- c("autm", "cube")
        sq$autm <- as.numeric(sq$autm)

        gr$autm[idx] <- sq$autm
        gr$cube[idx] <- sq$cube
        idx <- which(gr$cube == cube[2])
        strands[idx] <- "-"
        strand(gr) <- strands
    }
    return(gr)
}

digit_rep <- function(x,
    base = 2,
    ndigits) {
    if (any(x < 0)) {
        stop("*** 'x' must be non-negative integers")
    }
    if (!is.numeric(x)) {
        stop("*** 'x' must be integer-valued")
    } else {
        if (x != as.integer(x)) {
            stop("*** 'x' must be integer-valued")
        }
        x <- as.integer(x)
    }

    dig <- matrix(0, nrow = ndigits, ncol = length(x))
    if (ndigits >= 1) {
        for (i in ndigits:1) {
            dig[i, ] <- x %% base
            if (i > 1) {
                x <- x %/% base
            }
        }
    }
    return(dig)
}
