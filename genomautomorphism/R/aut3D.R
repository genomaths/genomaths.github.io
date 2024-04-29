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

#' @rdname aut3D
#' @aliases aut3D
#' @title Compute the Automorphisms of Mutational Events Between two Codon
#' Sequences Represented in Z5^3.
#' @description Given two codon sequences represented in the Z5^3 Abelian
#' group, this function computes the automorphisms describing codon mutational
#' events.
#' @details Automorphisms in Z5^3' are described as functions
#' \eqn{f(x) = A x mod Z5}, where A is diagonal matrix, as noticed in
#' reference (4).
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
#' @param field  A character string denoting the Galois field where the 3D
#' automorphisms are estimated. This can be 'GF(4)' or 'GF(5)', but only
#' 'GF(5)' is implemented so far.
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
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam multicoreWorkers
#' @importFrom methods new
#' @import Biostrings
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @references
#' \enumerate{
#'  \item Sanchez R, Morgado E, Grau R. Gene algebra from a genetic code
#'  algebraic structure. J Math Biol. 2005 Oct;51(4):431-57.
#'  doi: 10.1007/s00285-005-0332-8. Epub 2005 Jul 13. PMID: 16012800. (
#'  [PDF](https://arxiv.org/pdf/q-bio/0412033.pdf)).
#'  \item Robersy Sanchez, Jesus Barreto (2021) Genomic Abelian Finite
#'   Groups. \url{https://doi.org/10.1101/2021.06.01.446543}.
#'  \item M. V Jose, E.R. Morgado, R. Sanchez, T. Govezensky, The 24 possible
#'  algebraic representations of the standard genetic code in six or in three
#'  dimensions, Adv. Stud. Biol. 4 (2012) 119-152.[PDF](https://is.gd/na9eap).
#'  \item R. Sanchez. Symmetric Group of the Genetic-Code Cubes. Effect of the
#'  Genetic-Code Architecture on the Evolutionary Process MATCH Commun. Math.
#'  Comput. Chem. 79 (2018) 527-560. [PDF](https://bit.ly/2Z9mjM7).
#' }
#' @examples
#' ## Load a pairwise alignment
#' data("aln", package = "GenomAutomorphism")
#' aln
#'
#' ## Automorphism on Z5^3
#' autms <- aut3D(seq = aln)
#' autms
#'
aut3D <- function(seq = NULL,
    filepath = NULL,
    cube = c("ACGT", "TGCA"),
    cube_alt = c("CATG", "GTAC"),
    field = "GF5",
    start = NA,
    end = NA,
    chr = 1L,
    strand = "+",
    genetic_code = getGeneticCode("1"),
    num.cores = multicoreWorkers(),
    tasks = 0L,
    verbose = TRUE) {
    if (is.null(filepath) && is.null(seq)) {
        stop("*** One of the arguments 'seq' or 'filepath' must be given.")
    }

    if (!is.null(seq)) {
        if (!inherits(seq, c("DNAStringSet", "DNAMultipleAlignment"))) {
            stop(
                "*** Agument 'seq' must belong to 'DNAStringSet'",
                " DNAMultipleAlignment class."
            )
        }
        if (any(nchar(seq) %% 3 != 0)) {
            stop(
                "*** The argument of 'seq' must be a pairwise alignment",
                " of codon sequences."
            )
        }
    }

    autm1 <- automorfismos_3D(
        seq = seq,
        filepath = filepath,
        cube = cube,
        start = start,
        end = end,
        chr = chr,
        strand = strand,
        field = field,
        genetic_code = genetic_code,
        num.cores = num.cores,
        tasks = tasks,
        verbose = verbose
    )

    idx <- which(is.na(autm1$autm))

    if (length(idx) > 0) {
        autm2 <- automorfismos_3D(
            seq = seq,
            filepath = filepath,
            cube = cube_alt,
            start = start,
            end = end,
            chr = chr,
            strand = strand,
            field = field,
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



## ===================== Auxiliary function ===========================

automorfismos_3D <- function(seq,
    filepath,
    cube,
    start = NA,
    end = NA,
    chr = 1L,
    strand = "+",
    field,
    genetic_code = getGeneticCode("1"),
    num.cores,
    tasks,
    verbose) {
    seq <- get_coord(
        x = seq,
        output = "all",
        base_seq = FALSE,
        filepath = filepath,
        cube = cube[1],
        group = "Z5^3",
        start = start,
        end = end,
        chr = chr,
        strand = strand
    )

    gr <- seq@SeqRanges
    gr$aa1 <- translation(gr$seq1, genetic.code = genetic_code)
    gr$aa2 <- translation(gr$seq2, genetic.code = genetic_code)
    gr$coord1 <- seq@CoordList$coord1
    gr$coord2 <- seq@CoordList$coord2
    gr$autm <- "1,1,1"
    gr$cube <- cube[1]

    idx <- apply(seq@CoordList$coord1 != seq@CoordList$coord2, 1, any)
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

    ### -------------------------------------------------------------- #

    if (length(idx) != 0) {
        seq <- bplapply(idx, function(k) {
            c1 <- seq@CoordList$coord1[k, ]
            c2 <- seq@CoordList$coord2[k, ]

            s <- try(mapply(modeq, a = c1, b = c2, 
                            MoreArgs = list(n = 5L)),
                silent = TRUE
            )

            if (any(s == -1) || inherits(s, "try-error")) {
                s <- try(mapply(modeq, a = (5 - c1), b = (5 - c2), 
                                MoreArgs = list(n = 5L)),
                    silent = TRUE
                )
                if (!(any(s == -1) || inherits(s, "try-error"))) {
                    s <- c(paste0(s, collapse = ","), cube[2])
                }
            } else {
                s <- c(paste0(s, collapse = ","), cube[1])
            }
            if (any(s == -1) || inherits(s, "try-error")) {
                s <- c(0, "Trnl")
            }
            return(s)
        }, BPPARAM = bpparam)

        seq <- do.call(rbind, seq)
        seq <- data.frame(seq)
        colnames(seq) <- c("autm", "cube")
        gr$autm[idx] <- seq$autm
        gr$cube[idx] <- seq$cube
    }
    return(gr)
}
