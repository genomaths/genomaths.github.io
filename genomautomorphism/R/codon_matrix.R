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
##

## ==================== codon_matrix def. ===========================

#' Codon Coordinate Matrix
#' @rdname codon_matrix
#' @aliases codon_matrix
#' @import Biostrings
#' @description
#' This function build the coordinate matrix for each sequence from an aligned
#' set of DNA codon sequences.
#' @details
#' The purpose of this function is making the codon coordinates from multiple
#' sequence alignments (MSA) available for further downstream statistical
#' analyses, like those reported in references (1) and (2).
#' 
#' @param base A \code{\link[Biostrings]{DNAMultipleAlignment}}, a 
#' \code{\link[Biostrings]{DNAStringSet}}, or a [BaseSeqMatrix].
#' @author Robersy Sanchez <https://genomaths.com>
#' @returns A [ListCodonMatrix] class object with the codon coordinate on its
#' metacolumns.
#' @seealso [codon_coord], [base_coord] and [base2int].
#' @export
#' @references
#' \enumerate{
#'  \item Lorenzo-Ginori, Juan V., Aníbal Rodríguez-Fuentes, Ricardo Grau 
#'  Ábalo, and Robersy Sánchez Rodríguez. "Digital signal processing in the
#'  analysis of genomic sequences." Current Bioinformatics 4, no. 1 (2009):
#'  28-40.
#'  \item Sanchez, Robersy. "Evolutionary analysis of DNA-protein-coding
#'  regions based on a genetic code cube metric." Current Topics in Medicinal
#'  Chemistry 14, no. 3 (2014): 407-417.
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

#' @references 
#' 1. 
#' 2. 
#' @examples
#' ## Load the MSA of Primate BRCA1 DNA repair genes
#' data("brca1_aln")
#' 
#' ## Get the DNAStringSet for the first 33 codons and apply 'codon_matrix'
#' brca1 <- unmasked(brca1_aln)
#' brca1 <- subseq(brca1, start = 1, end = 33)
#' codon_matrix(brca1)
#' 
#' ## Get back the alignment object and apply 'codon_matrix' gives us the 
#' ## same result.
#' brca1 <- DNAMultipleAlignment(as.character(brca1))
#' codon_matrix(brca1)
#' 
setGeneric(
    "codon_matrix",
    function(base, ...) {
        standardGeneric("codon_matrix")
    }
)

## ==================== BaseSeqMatrix ========================

#' @aliases codon_matrix
#' @rdname codon_matrix
#' @param num.cores,tasks Parameters for parallel computation using package
#' \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#' use, i.e. at most how many child processes will be run simultaneously (see
#' \code{\link[BiocParallel]{bplapply}} and the number of tasks per job (only
#' for Linux OS).
#' @param verbose If TRUE, prints the function log to stdout
#' @param ... Not in use yet.
#' @import Biostrings
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam multicoreWorkers
#' @export
setMethod(
    "codon_matrix", signature(base = "BaseSeqMatrix"),
    function(
        base, 
        num.cores = 1L, 
        tasks = 0L,
        verbose = TRUE,
        ...) {
    
        group <- base@group
        cube <- base@cube
        alias <- base@seq_alias
        
        if (length(base) %% 3 != 0) {
            stop(
                "*** 'base' argument does not carry base-triplet sequences.",
                " A base-triplet sequence is multiple of 3."
            )
        }
        
        chr <- unique(as.character(seqnames(base)))
        strands <- unique(as.character(strand(base)))
        
        base <- mcols(base)
        sqnms <- colnames(base)
        idx_coord <- grep("S", sqnms)
        
        base <- data.frame(base[, idx_coord])
        colnames(base) <- sqnms
            
        f <- factor(as.vector(slapply(seq_len(nrow(base)/3), rep, times = 3)))
        base <- split(base, f)
        nms <- paste0("codon.", names(base))
        
        if (num.cores > 1L) {
        ## -------------- Setting parallel computation --------------- #
            dc <- multicoreWorkers()
            if (num.cores > dc)
                num.cores <- dc
            
            progressbar <- FALSE
            if (verbose)
                progressbar <- TRUE
            if (Sys.info()["sysname"] == "Linux")
                bpparam <- MulticoreParam(workers = num.cores, tasks = tasks,
                                        progressbar = progressbar)
            else
                bpparam <- SnowParam(workers = num.cores, type = "SOCK",
                                    progressbar = progressbar)
            ## ----------------------------------------------------------- #
            
            base <- bplapply(
                seq_along(base), 
                function(k) {
                    pos <- as.numeric(rownames(base[[k]]))
                    x <- data.frame(
                            seqnames = chr,
                            start = pos[1],
                            end = pos[3],
                            strand = strands,
                            base[[k]]
                        )
                    x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
                    CodonMatrix(x, cube = cube, group = group, 
                                seq_alias = alias[k])
                }, BPPARAM = bpparam
            )
        }
        else {
            for (k in seq_along(base)) {
                pos <- as.numeric(rownames(base[[k]]))
                x <- data.frame(
                    seqnames = chr,
                    start = pos[1],
                    end = pos[3],
                    strand = strands,
                    base[[k]]
                )
                x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
                base[[k]] <- CodonMatrix(x, cube = cube, group = group, 
                                        seq_alias = alias[k])
            }
        }
        names(base) <- nms
        
        base <- ListCodonMatrix(object = base, cube = cube, group = group,
                        seq_alias = alias, names = nms)
        
        return(base)
    }
)

## ==================== DNAStringSet ========================

#' @rdname codon_matrix
#' @aliases codon_matrix
#' @export
setMethod(
    "codon_matrix", signature(base = "DNAStringSet"),
    function(
        base, 
        cube = c("ACGT", "AGCT", "TCGA", "TGCA", "CATG", "GTAC", 
                "CTAG", "GATC", "ACTG", "ATCG", "GTCA", "GCTA", 
                "CAGT", "TAGC", "TGAC", "CGAT", "AGTC", "ATGC", 
                "CGTA", "CTGA", "GACT", "GCAT", "TACG", "TCAG"),
        group = c("Z4", "Z5"),
        num.cores = 1L, 
        tasks = 0L,
        verbose = TRUE) {

        base <- base_matrix(base = base, cube = cube, group = group)
        
        codon_matrix(base = base, num.cores = num.cores, tasks = tasks)
    }
)



## ==================== DNAMultipleAlignment ========================

#' @rdname codon_matrix
#' @aliases codon_matrix
#' @param cube A character string denoting one of the 24 Genetic-code cubes,
#' as given in references (3-4).
#' @param group A character string denoting the group representation for the
#' given base or codon as shown in reference (3-4).
#' @export
setMethod(
    "codon_matrix", signature(base = "DNAMultipleAlignment"),
    function(
        base, 
        cube = c("ACGT", "AGCT", "TCGA", "TGCA", "CATG", "GTAC", 
                "CTAG", "GATC", "ACTG", "ATCG", "GTCA", "GCTA", 
                "CAGT", "TAGC", "TGAC", "CGAT", "AGTC", "ATGC", 
                "CGTA", "CTGA", "GACT", "GCAT", "TACG", "TCAG"),
        group = c("Z4", "Z5"),
        num.cores = 1L, 
        tasks = 0L,
        verbose = TRUE) {
        
        base <- base@unmasked
        base <- base_matrix(base = base, cube = cube, group = group)
        
        codon_matrix(base = base, num.cores = num.cores, tasks = tasks)
    }
)

