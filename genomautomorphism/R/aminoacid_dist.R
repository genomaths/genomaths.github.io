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

#' Distance Between Aminoacids in Terms of Codon Distance
#' @rdname aminoacid_dist
#' @aliases aminoacid_dist
#' @description This function computes the distance between aminoacids in 
#' terms of a statistic of the corresponding codons. The possible statistics
#' are: 'mean', 'median', or some user defined function.
#' @details Only aminoacids sequences given in the following alphabet are 
#' accepted: "A","R","N","D","C","Q","E","G","H","I","L","K", "M","F","P",
#' "S","T","W","Y","V", "*", "-", and "X"; where symbols "*" and "-" denote 
#' the presence a stop codon and of a gap, respectively, and letter "X" 
#' missing information, which are then taken as a gap.
#' 
#' The distance between any aminoacid and any of the non-aminoacid symbols is
#' the ceiling of the greater distance found in the corresponding aminoacid
#' distance matrix. 
#' 
#' @param aa1,aa2 A character string of codon sequences, i.e., sequences of 
#' DNA base-triplets. If only 'x' argument is given, then it must be a
#' \code{\link[Biostrings]{DNAStringSet-class}} object.
#' @param weight A numerical vector of  weights to compute weighted Manhattan 
#' distance between codons. If \eqn{weight = NULL}, then 
#' \eqn{weight = (1/4,1,1/16)} for \eqn{group = "Z4"} and 
#' \eqn{weight = (1/5,1,1/25)} for \eqn{group = "Z5"} (see 
#' \code{\link{codon_dist}}). 
#' @param stat The name of some statistical function summarizing data like
#' 'mean', 'median', or some user defined function ('user_def'). If 
#' \eqn{stat = 'user_def'}, then function must have a logical argument named 
#' 'na.rm' addressed to remove missing (NA) data (see e.g., 
#' \code{\link[base]{mean}}).
#' @param genetic_code A single string that uniquely identifies the genetic 
#' code to extract. Should be one of the values in the id or name2 columns of 
#' \code{\link[Biostrings]{GENETIC_CODE_TABLE}}.
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
#' @importFrom utils data
#' @export
#' @seealso \code{\link{automorphisms}} and \code{\link{codon_coord}}
#' @references
#' \enumerate{
#'  \item Sanchez R. Evolutionary Analysis of DNA-Protein-Coding Regions Based
#'  on a Genetic Code Cube Metric. Curr. Top. Med. Chem. 2014;14: 407–417. 
#'  \url{https://doi.org/10.2174/1568026613666131204110022}.
#'  \item M. V Jose, E.R. Morgado, R. Sanchez, T. Govezensky, The 24 possible
#'  algebraic representations of the standard genetic code in six or in three
#'  dimensions, Adv. Stud. Biol. 4 (2012) 119-152.[PDF](https://is.gd/na9eap).
#'  \item R. Sanchez. Symmetric Group of the Genetic-Code Cubes. Effect of the
#'  Genetic-Code Architecture on the Evolutionary Process MATCH Commun. Math.
#'  Comput. Chem. 79 (2018) 527-560. [PDF](https://is.gd/ZY1Gx8).
#' }
#' @seealso \code{\link{codon_dist}}
#' @examples
#' ## Write down to aminoacid sequences
#' x <- "A*LTHMC"
#' y <- "AAMTDM-"
#' 
#' aminoacid_dist(aa1 = x, aa2 = y)
#' 
#' ## Let's create an AAStringSet-class object
#' aa <- AAStringSet(c(x, y))
#' 
#' aminoacid_dist(aa1 = aa)
#' 
#' ## Let's select cube "GCAT" and group "Z5"
#' aminoacid_dist(aa1 = aa, group = "Z5", cube = "TCGA")
#' 
setGeneric(
    "aminoacid_dist",
    function(
        aa1,
        aa2,
        ...) {
        standardGeneric("aminoacid_dist")
    }
)

## =============== Characters =====================


#' @rdname aminoacid_dist
#' @aliases aminoacid_dist
#' @import S4Vectors
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
setMethod(
    "aminoacid_dist", signature(aa1 = "character", aa2 = "character"),
    function(
        aa1,
        aa2,
        weight = NULL,
        stat = c("mean", "median", "user_def"),
        genetic_code = "1",
        group = c("Z4", "Z5"),
        cube = c("ACGT", "AGCT", "TCGA", "TGCA", "CATG", 
                "GTAC", "CTAG", "GATC", "ACTG", "ATCG", 
                "GTCA", "GCTA", "CAGT", "TAGC", "TGAC", 
                "CGAT", "AGTC", "ATGC", "CGTA", "CTGA", 
                "GACT", "GCAT", "TACG", "TCAG"),
        num.cores = 1L,
        tasks = 0L,
        verbose = FALSE) {
        
        cdm_z64 <- NULL
        group <- match.arg(group)
        cube <- match.arg(cube)
        stat <- match.arg(stat)
        
        if (genetic_code == "1" && group == "Z4") {
            data("cdm_z64", package = "GenomAutomorphism", 
                envir = environment(), overwrite = TRUE)
            cdm <- cdm_z64[[ cube ]]
        }
        else
            cdm <- codon_dist_matrix(genetic_code = genetic_code, 
                                    cube = cube, weight = weight, 
                                    group = group, output = "vector",
                                    num.cores = num.cores)
        
        alf <- c("A","R","N","D","C","Q","E","G","H","I","L","K",
                "M","F","P","S","T","W","Y","V", "*", "-", "X")
        
        stat <- switch(stat,
                    "mean" = mean,
                    "median" = median
        )
        
        if (!is.function(stat))
            stop("*** Argument 'stat' must a function.")
        
        if (length(aa1) == 1) {
            if (nchar(aa1) > 1) 
                aa1 <- str2chr(aa1)
        }

        if (length(aa2) == 1) {
            if (nchar(aa2) > 1)
                aa2 <- str2chr(aa2)
        }
        
        if (any(!is.element(unique(aa1), alf)))
            stop("*** Argument 'aa1' carries letters outside the",
                " aminoacid alphabet.")
        if (any(!is.element(unique(aa2), alf)))
            stop("*** Argument 'aa2' carries letters outside the",
                " aminoacid alphabet.")
        
        gc <- getGeneticCode(id_or_name2 = genetic_code)
        nms <- names(gc)
    
        max_dist <- ceiling(max(cdm))
        
        if (length(aa1) == 1) {
            if (any(is.element(c(aa1, aa2), c("X", "-", "*")))) {
                dm <- max_dist
            }
            else {
                aa1 <- match(aa1, gc)
                aa2 <- match(aa2, gc)
                
                cd1.cd2 <- c(as.vector(outer(nms[aa1], nms[aa2],
                                            FUN = paste, sep = ".")),
                            as.vector(outer(nms[aa2], nms[aa1], 
                                            FUN = paste, sep = ".")))

                dm <- try(stat(cdm[ na.omit(match(cd1.cd2, names(cdm))) ],
                            na.rm = TRUE),
                        silent = TRUE)
                if (inherits(dm, "try-error"))
                    stop("*** The statistic given in function 'stat'",
                        " cannot be computes.")
            }
        }
        else {
            # ## ---------- Setting parallel computation ------------ #
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
            # ## ---------------------------------------------------- #
            
            dm <- bplapply(seq_along(aa1), function(k) {
                
                if (any(is.element(c(aa1[k], aa2[k]), c("X", "-", "*")))) {
                    dm <- max_dist
                }
                else {
                    a1 <- grep(aa1[k], gc)
                    a2 <- grep(aa2[k], gc)
                    
                    if (length(nms[a1]) > 1 || length(nms[a2]) > 1) {
                        cd1.cd2 <- c(as.vector(outer(nms[a1], nms[a2],
                                                    FUN = paste, sep = ".")),
                                    as.vector(outer(nms[a2], nms[a1], 
                                                    FUN = paste, sep = ".")))
                    }
                    else {
                        if (nms[a1] != nms[a2])
                            cd1.cd2 <- paste(nms[a1], nms[a2], sep = ".")
                        else
                            cd1.cd2 <- NA
                    }
                    
                    if (all(!is.na(cd1.cd2))) {
                        dm <- try(stat(
                                cdm[na.omit(match(cd1.cd2, names(cdm)))],
                                    na.rm = TRUE),
                                silent = TRUE)
                    }
                    else
                        dm <- 0
                    if (inherits(dm, "try-error"))
                        stop("*** The statistic given in function 'stat'",
                            " cannot be computes.")
                }
                return(dm)
            })
            names(dm) <- paste(aa1, aa2, sep = ".")
            dm <- unlist(dm)
        }
        return(dm)
    }
)

## =============== DNAStringSet =====================

#' @rdname aminoacid_dist
#' @aliases aminoacid_dist
#' @import S4Vectors
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
setMethod(
    "aminoacid_dist", signature(aa1 = "DNAStringSet"),
    function(
        aa1,
        weight = NULL,
        stat = c("mean", "median", "user_def"),
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
        stat <- match.arg(stat)
        
        aa1 <- translation(x = aa1)
        aa1 <- as.character(aa1)

        aa1 <- aminoacid_dist(aa1 = aa1[1], aa2 = aa1[2], stat = stat, 
                            weight = weight, group = group, cube = cube, 
                            num.cores = num.cores, tasks = tasks, 
                            verbose = verbose)
        return(aa1)
    }
)


## =============== AAStringSet =====================


#' @rdname aminoacid_dist
#' @aliases aminoacid_dist
#' @import S4Vectors
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
setMethod(
    "aminoacid_dist", signature(aa1 = "AAStringSet"),
    function(
        aa1,
        weight = NULL,
        stat = c("mean", "median", "user_def"),
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
        stat <- match.arg(stat)
        
        aa1 <- as.character(aa1)
        
        aa1 <- aminoacid_dist(aa1 = aa1[1], aa2 = aa1[2], 
                    weight = weight, stat = stat, group = group, 
                    cube = cube, num.cores = num.cores, 
                    tasks = tasks,  verbose = verbose)
        return(aa1)
    }
)


## =============== CodonGroup_OR_Automorphisms =====================

#' @rdname aminoacid_dist
#' @aliases aminoacid_dist
#' @import S4Vectors
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
setMethod(
    "aminoacid_dist", signature(aa1 = "CodonGroup_OR_Automorphisms"),
    function(
        aa1,
        weight = NULL,
        stat = c("mean", "median", "user_def"),
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
        stat <- match.arg(stat)
        
        aa1 <- aminoacid_dist(aa1 = aa1$aa1, aa2 = aa1$aa2,
                    weight = weight, stat = stat, group = group, 
                    cube = cube, num.cores = num.cores,
                    tasks = tasks, verbose = verbose)
        return(aa1)
    }
)

