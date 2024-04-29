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

#' @rdname automorphisms
#' @title Compute the Automorphisms of Mutational Events Between two Codon
#' Sequences Represented in a Given Abelian group.
#' @description Given two codon sequences represented in a given Abelian
#' group, this function computes the automorphisms describing codon mutational
#' events. Basically, this function is a wrapping to call the corresponding
#' function for a specified Abelian group.
#'
#' @details Herein, automorphisms are algebraic descriptions of mutational
#' event observed in codon sequences represented on different Abelian groups.
#' In particular, as described in references (3-4), for each representation of
#' the codon set on a defined Abelian group there are 24 possible isomorphic
#' Abelian groups. These Abelian groups can be labeled based on the DNA
#' base-order used to generate them. The set of 24 Abelian groups can be
#' described as a group isomorphic to the symmetric group of degree four
#' (\eqn{S_4}, see reference (4)). Function \code{\link{automorphismByRanges}}
#' permits the classification of the pairwise alignment of protein-coding
#' sub-regions based on the mutational events observed on it and on the
#' genetic-code cubes that describe them.
#'
#' Automorphisms in Z5, Z64 and Z125 are described as functions
#' \eqn{f(x) = k x mod 64} and \eqn{f(x) = k x mod 125}, where k and x are
#' elements from the set of integers modulo 64 or modulo 125, respectively. If
#' an automorphisms cannot be found on any of the cubes provided in the
#' argument \eqn{cube}, then function \code{\link{automorphisms}} will search
#' for automorphisms in the cubes provided in the argument \eqn{cube_alt}.
#'
#' Automorphisms in Z5^3' are described as functions \eqn{f(x) = Ax mod Z5},
#' where A is diagonal matrix.
#'
#' Arguments \emph{\strong{cube}} and \emph{\strong{cube_alt}} must be
#' pairs of' dual cubes (see section 2.4 from reference 4).
#'
#' @param seqs An object from a \code{\link[Biostrings]{DNAStringSet}} or
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
#' call \eqn{dual cubes} and, as shown in reference (3), each pair integrates
#' group.
#' @param nms Optional. Only used if the DNA sequence alignment provided
#' carries more than two sequences. A character string giving short names for
#' the alignments to be compared. If not given then the automorphisms between
#' pairwise alignment are named as: 'aln_1', 'aln_2', and so on.
#' @param start,end,chr,strand Optional parameters required to build a
#' \code{\link[GenomicRanges]{GRanges-class}}. If not provided the default
#' values given for the function definition will be used.
#' @param group A character string denoting the group representation for the
#' given base or codon as shown in reference (1).
#' @param num.cores,tasks Parameters for parallel computation using package
#' \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#' use, i.e. at most how many child processes will be run simultaneously (see
#' \code{\link[BiocParallel]{bplapply}} and the number of tasks per job (only
#' for Linux OS).
#' @param verbose If TRUE, prints the progress bar.
#' @param ... Not in use.
#' @return This function returns a \code{\link{Automorphism-class}} object
#' with four columns on its metacolumn named: \emph{seq1}, \emph{seq2},
#' \emph{autm}, and \emph{cube}.
#'
#' @section Methods:
#' ### \code{\link{automorphismByRanges}}:
#' This function returns a \code{\link[GenomicRanges]{GRanges-class}} object.
#' Consecutive mutational events (on the codon sequence) described by
#' automorphisms on a same cube are grouped in a range.
#'
#' ### \code{\link{automorphism_bycoef}}
#' This function returns a \code{\link[GenomicRanges]{GRanges-class}} object.
#' Consecutive mutational events (on the codon sequence) described by
#' the same automorphisms coefficients are grouped in a range.
#'
#' ### \code{\link{getAutomorphisms}}
#' This function returns an AutomorphismList-class object as a list of
#' Automorphism-class objects, which inherits from
#' \code{\link[GenomicRanges]{GRanges-class}} objects.
#'
#' ### \code{\link{conserved_regions}}
#' Returns a \code{\link{AutomorphismByCoef}} class object containing the
#' requested regions.
#'
#' @seealso \code{\link{autZ64}}.
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
#' data("aln", package = "GenomAutomorphism")
#' aln
#'
#' ## Automorphism on "Z5^3"
#' autms <- automorphisms(seqs = aln, group = "Z5^3", verbose = FALSE)
#' autms
#'
#' ## Automorphism on "Z64"
#' autms <- automorphisms(seqs = aln, group = "Z64", verbose = FALSE)
#' autms
#'
#' ## Automorphism on "Z64" from position 1 to 33
#' autms <- automorphisms(
#'     seqs = aln,
#'     group = "Z64",
#'     start = 1,
#'     end = 33,
#'     verbose = FALSE
#' )
#' autms
#'
#' @aliases automorphisms
#' @export
setGeneric(
    "automorphisms",
    function(seqs = NULL,
    filepath = NULL,
    group = "Z4",
    ...) {
        standardGeneric("automorphisms")
    }
)

#' @aliases automorphisms
#' @rdname automorphisms
#' @importFrom foreach foreach %dopar% %:%
#' @importFrom doParallel registerDoParallel
#' @importFrom stats setNames
#' @importFrom parallel makeCluster stopCluster
#' @importFrom BiocParallel multicoreWorkers
#' @importFrom BiocGenerics width
#' @import Biostrings
#' @import S4Vectors
#' @export
setMethod(
    "automorphisms", signature(seqs = "DNAStringSet_OR_NULL"),
    function(seqs = NULL,
    filepath = NULL,
    group = c("Z5", "Z64", "Z125", "Z5^3"),
    cube = c("ACGT", "TGCA"),
    cube_alt = c("CATG", "GTAC"),
    nms = NULL,
    start = NA,
    end = NA,
    chr = 1L,
    strand = "+",
    num.cores = multicoreWorkers(),
    tasks = 0L,
    verbose = TRUE) {
        if (is.null(seqs) && is.null(filepath)) {
            stop("*** One of the arguments 'seqs' or 'filepath' must
                be provided.")
        }

        group <- match.arg(group)
        nr <- NA

        if (!is.null(filepath) && is.character(filepath)) {
            seqs <- readDNAMultipleAlignment(filepath = filepath)
        }

        if (inherits(seqs, "DNAStringSet")) {
            nr <- length(seqs)
        }

        if (inherits(seqs, "DNAMultipleAlignment")) {
            nr <- nrow(seqs)
        }

        if (is.na(nr)) {
            stop("*** Proper argument 'filepath' or 'seqs' must be provided.")
        }

        if (nr < 3) {
            seqs <- selectAutomorphism(
                seq = seqs,
                filepath = NULL,
                group = group,
                cube = cube,
                cube_alt = cube_alt,
                start = start,
                end = end,
                chr = chr,
                strand = strand,
                num.cores = num.cores,
                tasks = tasks,
                verbose = verbose
            )
        } else {
            ## just to get the ranges
            gr <- seqs@unmasked[c(1, 1)]
            gr <- selectAutomorphism(
                seq = gr,
                filepath = NULL,
                group = group,
                cube = cube,
                cube_alt = cube_alt,
                start = start,
                end = end,
                chr = chr,
                strand = strand,
                num.cores = num.cores,
                tasks = tasks,
                verbose = verbose
            )
            mcols(gr) <- NULL

            ## ------------ Setting up parallel computation ------------ #

            num.cores <- floor(num.cores / 2)
            no_cores <- num.cores
            
            if (Sys.info()["sysname"] == "Linux") 
                cl <- makeCluster(num.cores, type = "FORK")
            else
                cl <- makeCluster(num.cores, type = "SOCK")
            
            registerDoParallel(cl)

            if (is.null(nms)) {
                nms <- paste0("aln_", seq_len(nr))
            }
            nams <- function(x) {
                nms <- rev(nms[(length(nms) - seq_along(x) + 1)])
                nms
            }

            ## -------------------------------------------------------- #
            k <- j <- NULL
            seqs <- try(foreach(
                k = seq_len(nr - 1),
                .final = function(x) setNames(x, nms[-nr])
            ) %:%
                foreach(
                    j = seq((k + 1), nr, 1),
                    .final = function(x) setNames(x, nams(x))
                ) %dopar% {
                    aln <- seqs@unmasked[c(k, j)]

                    aln <- selectAutomorphism(
                        seq = aln,
                        filepath = NULL,
                        group = group,
                        cube = cube,
                        cube_alt = cube_alt,
                        start = start,
                        end = end,
                        chr = chr,
                        strand = strand,
                        num.cores = num.cores,
                        tasks = tasks,
                        verbose = verbose
                    )
                    mcols(aln)
                }, silent = TRUE)

            if (inherits(seqs, "try-error")) {
                stopCluster(cl)
                stop("*** Automorphism cannot be computed from
                    the MSA.")
            } else {
                stopCluster(cl)
                seqs <- unlist(seqs, recursive = FALSE)
                seqs <- as.AutomorphismList(seqs, grs = gr)
            }
        }
        return(seqs)
    }
)


## ===================== Auxiliary functions ===========================

selectAutomorphism <- function(seq,
    filepath,
    group,
    cube,
    cube_alt,
    start,
    end,
    chr,
    strand,
    num.cores,
    tasks,
    verbose) {
    seq <- switch(group,
        "Z5" = autZ5(
            seq = seq,
            filepath = filepath,
            cube = cube,
            cube_alt = cube_alt,
            start = start,
            end = end,
            chr = chr,
            strand = strand,
            num.cores = num.cores,
            tasks = tasks,
            verbose = verbose
        ),
        "Z64" = autZ64(
            seq = seq,
            filepath = filepath,
            cube = cube,
            cube_alt = cube_alt,
            start = start,
            end = end,
            chr = chr,
            strand = strand,
            num.cores = num.cores,
            tasks = tasks,
            verbose = verbose
        ),
        "Z5^3" = aut3D(
            seq = seq,
            filepath = filepath,
            cube = cube,
            cube_alt = cube_alt,
            start = start,
            end = end,
            chr = chr,
            strand = strand,
            num.cores = num.cores,
            tasks = tasks,
            verbose = verbose
        ),
        "Z125" = autZ125(
            seq = seq,
            filepath = filepath,
            cube = cube,
            cube_alt = cube_alt,
            start = start,
            end = end,
            chr = chr,
            strand = strand,
            num.cores = num.cores,
            tasks = tasks,
            verbose = verbose
        )
    )
    return(seq)
}
