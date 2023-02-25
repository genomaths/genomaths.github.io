## Copyright (C) 2022 Robersy Sanchez <https://genomaths.com/>
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

#' @rdname get_mutscore
#' @title Get Mutation Score from an AAindex Matrix
#' @description This function is applied to get the mutation or contact 
#' potential scores representing the similarity/distance between amino acids
#' corresponding to substitution mutations. The score are retrieve from a
#' mutation matrix or a statistical protein contact potentials matrix from
#' \href{https://www.genome.jp/aaindex/}{AAindex} (ver.9.2).
#' @param aa1,aa2 A simple character representing an amino acids or a 
#' character string of letter from the amino acid alphabet or base-triplets
#' from the DNA/RNA alphabet.
#' @param acc Accession id for a specified mutation or contact potential 
#' matrix.
#' @param aaindex Database where the requested accession id is locate. The 
#' possible values are:  "aaindex2" or "aaindex3".
#' @param mutmat A mutation or any score matrix provided by the user.
#' @param alphabet Whether the alphabet is from the 20 amino acid (AA) or
#' four (DNA)/RNA base alphabet. This would prevent mistakes, i.e., 
#' the strings "ACG" would be a base-triplet on the DNA alphabet or simply
#' the amino acid sequence of alanine, cysteine, and glutamic acid.
#' @param num.cores,tasks Parameters for parallel computation using package
#' \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#' use, i.e. at most how many child processes will be run simultaneously
#' (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#' (only for Linux OS).
#' @param verbose If TRUE, prints the function log to stdout.
#' @param ... Not in use.
#' 
#' @details If a score matrix is provided by the user, then it must be a
#' symmetric matrix 20x20.
#' @seealso \code{\link{aa_mutmat}}, \code{\link{aaindex2}} and 
#' \code{\link{aaindex3}}
#' @author Robersy Sanchez <https://genomaths.com>
#' @export
#' @return A single numeric score or a numerical vector.
#' @examples 
#' ## Load the mutation matrices from database from the packages
#' data("aaindex2", package = "GenomAutomorphism" )
#' 
#' ## A single amino acids substitution mutation
#' get_mutscore("A", "C", acc = "MIYS930101", aaindex = aaindex2)
#' 
#' ## A tri-peptide mutation
#' get_mutscore(aa1 = "ACG", aa2 = "ATG", acc = "MIYS930101",
#'             aaindex = aaindex2, alphabet = "AA")
#' 
#' ## A single base-triple mutation, i.e., a single amino acid substitution 
#' ## mutation
#' get_mutscore(aa1 = "ACG", aa2 = "CTA", acc = "MIYS930101",
#'             aaindex = aaindex2, alphabet = "DNA")
#' 
#' ## Peptides can be also written as:
#' get_mutscore(aa1 = c("A","C","G"), aa2 = c("C","T","A"), 
#'             acc = "MIYS930101", aaindex = aaindex2, alphabet = "AA")

setGeneric("get_mutscore",
    function(
        aa1, 
        aa2, 
        ...) standardGeneric("get_mutscore"))

#' @aliases get_mutscore
#' @rdname get_mutscore
#' @export
setMethod("get_mutscore", signature(aa1 = "character", aa2 = "character"),
    function(
        aa1, 
        aa2, 
        acc = NULL,
        aaindex = NULL,
        mutmat = NULL,
        alphabet = c("AA", "DNA"),
        num.cores = 1L,
        tasks = 0L,
        verbose = FALSE,
        ...) {
        
        alphabet <- match.arg(alphabet)
        
        aa <- c("A","R","N","D","C","Q","E","G","H",
                "I","L","K","M","F","P","S","T","W","Y","V")
        dna <- c("A","C","G","T","U")
        
        if (length(aa1) != length(aa2))
            stop("Arguments 'aa1' & 'aa2' must have the same length.")
        
        if (length(aa1) == 1 && alphabet == "DNA" && nchar(aa1) > 1 && 
            nchar(aa1) != 3) {
            aa1 <- base2codon(aa1)
            aa2 <- base2codon(aa2)
            aa1 <- translation(aa1)
            aa2 <- translation(aa2)
            alphabet <- "AA"
        }
                
        if (length(aa1) == 1 && nchar(aa1) == 3 && alphabet == "DNA") {
            query <- unique(unlist(str2chr(c(aa1, aa2))))
            is.triplet <- unique(nchar(c(aa1, aa2))) == 3
            
            if (is.triplet && all(is.element(query, dna))) {
                aa1 <- translation(aa1)
                aa2 <- translation(aa2)
                alphabet <- "AA"
            }
        }
        
        if (length(aa1) == 1 && nchar(aa1) > 1 && alphabet == "AA") {
            aa1 <- str2chr(aa1)
            aa2 <- str2chr(aa2)
            if (length(aa1) != length(aa2))
                stop("Arguments 'aa1' annd 'aa2' must have the same",
                    "the number of characters.")
        }
        
        if (length(aa1) == 1 && (nchar(aa1) == 1)) {
            if (nchar(aa1) != nchar(aa2))
                stop("Arguments 'aa1' annd 'aa2' must have the same",
                    "the number of characters.")
            if (!is.null(mutmat))
                return(mutmat[ match(aa1, aa), match(aa2, aa) ])
            else {
                if (is.null(acc))
                    stop("The accesion ID for the mutation",
                        " matrix must be provides")
                if (is.null(aaindex))
                    stop("The name of amino acid index database",
                        " matrix must be provides")
                mutmat <- aa_mutmat(acc = acc, aaindex = aaindex)
                return(mutmat[ match(aa1, aa), match(aa2, aa) ])
            }
        }
        else {
            ch1 <- unique(nchar(aa1))
            ch2 <- unique(nchar(aa2))
            
            if (length(ch1) > 1)
                stop("All the element from 'aa1' must have the same number", 
                    " of characters.")
            if (any(ch1 != ch2)) 
                stop("Not all the elements of 'aa1' and 'aa2' has the same ",
                    "number of characters.")
            if (any(!is.element(ch1, c(1,3))))
                stop("The elements from 'aa1' must be single letters of ",
                    "amino acid laphabet or base-triples on DNA/RNA", 
                    " alphabets.")
                
            if (any(!is.element(ch2, c(1,3))))
                stop("The elements from 'aa2' must be single letters of ",
                    "amino acid laphabet or base-triples on DNA/RNA",
                    " alphabets.")
            if (num.cores == 1) {
                score <- slapply(seq_along(aa1), 
                        function(k) {
                            if (nchar(aa1[k]) != nchar(aa2[k]))
                                stop("Arguments aa1[k] & aa2[k] with ",
                                    "different number of characters.",
                                    " k=", k)
                            get_mutscore(
                                aa1 = aa1[k], 
                                aa2 = aa2[k], 
                                acc = acc,
                                aaindex = aaindex,
                                mutmat = mutmat,
                                alphabet = alphabet)
                        })
            }
            else {
                ## ------------ Setting parallel computing ------------- ##
                progressbar <- FALSE
                if (verbose) progressbar <- TRUE
                if (Sys.info()["sysname"] == "Linux") {
                    bpparam <- MulticoreParam(
                        workers = num.cores,
                        tasks = tasks,
                        progressbar = progressbar)
                } else {
                    bpparam <- SnowParam(
                        workers = num.cores,
                        type = "SOCK",
                        progressbar = progressbar)
                }
                ## ----------------------------------------------------- ##
                
                score <- bplapply(
                        seq_along(aa1), 
                            function(k) {
                                get_mutscore(
                                    aa1 = aa1[k], 
                                    aa2 = aa2[k], 
                                    acc = acc,
                                    aaindex = aaindex,
                                    mutmat = mutmat,
                                    alphabet = alphabet) 
                            },
                            BPPARAM = bpparam)
            }
            return(unlist(score))
        }
    }
)
