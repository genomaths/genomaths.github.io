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


#' Weighted Manhattan Distance Between Codons
#' @rdname codon_dist
#' @aliases codon_dist
#' @description This function computes the weighted Manhattan distance between
#' codons from two sequences as given in reference (1). That is, given two
#' codons \eqn{x} and \eqn{y} with coordinates on the set of integers modulo 5
#' ("Z5"): \eqn{x = (x_1, x_2, x_3)} and  \eqn{x = (y_1, y_2, y_3)} (see (1)),
#' the Weighted Manhattan distance between this two codons is defined as: 
#' 
#' \deqn{d_w(x,y) = |x_1 - y_1|/5 + |x_2 - y_2| + |x_3 -y_3|/25}
#' 
#' If the codon coordinates are given on "Z4", then the Weighted Manhattan
#' distance is define as:
#' 
#' \deqn{d_w(x,y) = |x_1 - y_1|/4 + |x_2 - y_2| + |x_3 -y_3|/16}
#' 
#' Herein, we move to the generalized version given in reference (3), for 
#' which:
#' 
#' \deqn{d_w(x,y) = |x_1 - y_1| w_1  + |x_2 - y_2| w_2 + |x_3 -y_3| w_3}
#' 
#' where we use the vector of \eqn{weight = (w_1, w_2, w_3)}.
#' 
#' @param x,y  A character string of codon sequences, i.e., sequences of DNA 
#' base-triplets. If only 'x' argument is given, then it must be a
#' \code{\link[Biostrings]{DNAStringSet-class}} object.
#' @param weight A numerical vector of  weights to compute weighted Manhattan 
#' distance between codons. If \eqn{weight = NULL}, then 
#' \eqn{weight = (1/4,1,1/16)} for \eqn{group = "Z4"} and 
#' \eqn{weight = (1/5,1,1/25)} for \eqn{group = "Z5"}.
#' @param group A character string denoting the group representation for the
#' given codon sequence as shown in reference (2-3).
#' @param cube A character string denoting one of the 24 Genetic-code cubes,
#' as given in references (2-3).
#' @param num.cores,tasks Parameters for parallel computation using package
#' \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#' use, i.e. at most how many child processes will be run simultaneously (see
#' \code{\link[BiocParallel]{bplapply}} and the number of tasks per job (only
#' for Linux OS).
#' @param verbose If TRUE, prints the progress bar.
#' @param ... Not in use yet.
#' @return A numerical vector with the pairwise distances between codons in 
#' sequences 'x' and 'y'.
#' @export
#' @seealso \code{\link{codon_dist_matrix}}, \code{\link{automorphisms}}, 
#' \code{\link{codon_coord}}, and \code{\link{aminoacid_dist}}.
#' @references
#' \enumerate{
#'  \item Sanchez R. Evolutionary Analysis of DNA-Protein-Coding Regions Based
#'  on a Genetic Code Cube Metric. Curr Top Med Chem. 2014;14: 407–417. 
#'  \url{https://doi.org/10.2174/1568026613666131204110022}.
#'  \item M. V Jose, E.R. Morgado, R. Sanchez, T. Govezensky, The 24 possible
#'  algebraic representations of the standard genetic code in six or in three
#'  dimensions, Adv. Stud. Biol. 4 (2012) 119-152.[PDF](https://is.gd/na9eap).
#'  \item R. Sanchez. Symmetric Group of the Genetic-Code Cubes. Effect of the
#'  Genetic-Code Architecture on the Evolutionary Process MATCH Commun. Math.
#'  Comput. Chem. 79 (2018) 527-560. [PDF](https://is.gd/ZY1Gx8).
#' }
#' @examples
#' ## Let's write two small DNA sequences
#' x = "ACGCGTGTACCGTGACTG"
#' y = "TGCGCCCGTGACGCGTGA"
#' 
#' codon_dist(x, y, group = "Z5")
#' 
#' ## Alternatively, data can be vectors of codons, i.e., vectors of DNA 
#' ## base-triplets (including gaps simbol "-").
#' x = c("ACG","CGT","GTA","CCG","TGA","CTG","ACG")
#' y = c("TGC","GCC","CGT","GAC","---","TGA","A-G")
#' 
#' ## Gaps are not defined on "Z4"
#' codon_dist(x, y, group = "Z4")
#' 
#' ## Gaps are considered on "Z5"
#' codon_dist(x, y, group = "Z5")
#' 
#' ## Load an Automorphism-class object
#' data(autm, package = "GenomAutomorphism")
#' codon_dist(x = head(autm,20), group = "Z4")
#' 
#' ## Load a pairwise alignment
#' data(aln, package = "GenomAutomorphism")
#' aln
#' 
#' codon_dist(x = aln, group = "Z5")
#' 
setGeneric(
    "codon_dist",
    function(
        x,
        y,
        ...) {
        standardGeneric("codon_dist")
    }
)

#' @rdname codon_dist
#' @aliases codon_dist
#' @import S4Vectors
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
setMethod(
    "codon_dist", signature(x = "DNAStringSet"),
    function(
        x, 
        weight = NULL,
        group = c("Z4", "Z5"),
        cube = c("ACGT", "AGCT", "TCGA", "TGCA", "CATG", 
                "GTAC", "CTAG", "GATC", "ACTG", "ATCG", 
                "GTCA", "GCTA", "CAGT", "TAGC", "TGAC", 
                "CGAT", "AGTC", "ATGC", "CGTA", "CTGA", 
                "GACT", "GCAT", "TACG", "TCAG"),
        num.cores = 1L,
        tasks = 0L,
        verbose = FALSE) {
        
        group <- match.arg(group)
        cube <- match.arg(cube)
        ls <- unique(width(x))
        if (length(ls) > 1)
            stop("*** Codon sequences 'x' and 'y' must be aligned.")
        if ((ls %% 3) != 0)
            stop("*** Arguments 'x' and 'y' are not codon sequences.")
        
        if (is.null(weight)) {
            if (group == "Z4")
                weight <- c(1/4, 1, 1/16)
            else
                weight <- c(1/5, 1, 1/25)
        }
        
        x <- base_coord(base = x, cube = cube, group = group )
        x <- mcols(x)
        
        idx_coord <- grep("coord", colnames(x))
        
        x <- data.frame(x[, idx_coord])
        f <- factor(as.vector(slapply(seq_len(nrow(x) / 3), rep, times = 3)))
        x <- split(x, f)
    
        if (ls > 1 && num.cores == 1) {
            x <- slapply(seq(x), 
                    function(k) {
                        cds <- x[[k]]
                        weighted_manhattan(x = cds[, 1], y = cds[, 2],
                                        w = weight, group = group)
                    })
        }
        
        if (ls > 1 && num.cores > 1) {
            
            # ## -------------- Setting parallel computation ------------- #
            
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
            
            # ## -------------------------------------------------------- #
            
            x <- bplapply(seq(x), 
                    function(k) {
                        cds <- x[[k]]
                        weighted_manhattan(x = cds[, 1], y = cds[, 2],
                                    w = weight, group = group)
                    }, BPPARAM = bpparam)
            x <- unlist(x)
        }
        return(x)
    }
)

#' @rdname codon_dist
#' @aliases codon_dist
#' @import S4Vectors
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
setMethod(
    "codon_dist", signature(x = "character"),
    function(
        x, 
        y,
        weight = NULL,
        group = c("Z4", "Z5"),
        cube = c("ACGT", "AGCT", "TCGA", "TGCA", "CATG", 
                "GTAC", "CTAG", "GATC", "ACTG", "ATCG", 
                "GTCA", "GCTA", "CAGT", "TAGC", "TGAC", 
                "CGAT", "AGTC", "ATGC", "CGTA", "CTGA", 
                "GACT", "GCAT", "TACG", "TCAG"),
        num.cores = 1L,
        tasks = 0L,
        verbose = FALSE) {
        
        group <- match.arg(group)
        cube <- match.arg(cube)
        
        if (length(x) > 1) {
            x <- c(paste0(x, collapse = ""))
            y <- c(paste0(y, collapse = ""))
        }
        
        x <- DNAStringSet(c(x, y))
        x <- codon_dist(x = x, weight = weight, group = group, cube = cube, 
                num.cores = num.cores, tasks = tasks, verbose = verbose)
        return(x)
    }
)


setClassUnion("CodonGroup_OR_Automorphisms", 
            c("CodonGroup", "Automorphism", "AutomorphismByCoef"))

#' @rdname codon_dist
#' @aliases codon_dist
#' @import S4Vectors
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
setMethod(
    "codon_dist", signature(x = "CodonGroup_OR_Automorphisms"),
    function(
        x, 
        weight = NULL,
        group = c("Z4", "Z5"),
        cube = c("ACGT", "AGCT", "TCGA", "TGCA", "CATG", 
                "GTAC", "CTAG", "GATC", "ACTG", "ATCG", 
                "GTCA", "GCTA", "CAGT", "TAGC", "TGAC", 
                "CGAT", "AGTC", "ATGC", "CGTA", "CTGA", 
                "GACT", "GCAT", "TACG", "TCAG"),
        num.cores = 1L,
        tasks = 0L,
        verbose = FALSE) {
        
        group <- match.arg(group)
        cube <- match.arg(cube)
        
        x <- DNAStringSet(c(paste0(x$seq1, collapse = ""), 
                            paste0(x$seq2, collapse = "")))

        x <- codon_dist(x = x, weight = weight, group = group, cube = cube, 
                        num.cores = num.cores, tasks = tasks, 
                        verbose = verbose)
        return(x)
    }
)



## ------------------------- Auxiliary functions --------------------------

weighted_manhattan <- function(x, y, w, group) {
    if (group == "Z4") {
        dst <- abs(x[1] - y[1]) * w[1]  + abs(x[2] - y[2]) * w[3] + 
                    abs(x[3] - y[3]) * w[3]
    }else {
        dst <- abs(x[1] - y[1]) * w[1] + abs(x[2] - y[2]) * w[2] + 
                    abs(x[3] - y[3]) * w[3]
    }
    return(dst)
}

