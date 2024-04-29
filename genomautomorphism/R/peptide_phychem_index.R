## Copyright (C) 2022-2024 Robersy Sanchez <https://genomaths.com/>
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


## =========================== Definition ==========================

#' @rdname peptide_phychem_index
#' @aliases peptide_phychem_index
#' @title Amino acid numerical matrix
#' @description 
#' This function applies numerical indices representing various physicochemical
#' and biochemical properties of amino acids and pairs of amino acids to DNA
#' protein-coding or to aminoacid sequences. As results, DNA protein-coding or
#' the aminoacid sequences are represented as numerical vectors which can be
#' subject of further downstream statistical analysis and digital signal
#' processing.
#' 
#' @details
#' If a DNA sequence is given, then it is assumed that it is a DNA base-triplet
#' sequence, i.e., the base sequence must be multiple of 3. 
#' 
#' Errors can be originated if the given sequences carry letter which are not
#' from the DNA or aminoacid alphabet.
#' 
#' 
#' @examples
#' ## Let's create DNAStringSet-class object
#' base <- DNAStringSet(x = c( seq1 ='ACGTCATAAAGT',
#'                             seq2 = 'GTGTAATACAGT',
#'                             seq3 = 'TCCTCATAAGGT'))
#' 
#' ## The stop condon 'TAA' yields NA
#' aa <- peptide_phychem_index(base, acc = "EISD840101")
#' aa
#' 
#' ## Description of the physicochemical index
#' slot(aa, 'phychem')
#' 
#' ## Get the aminoacid sequences. The stop codon 'TAA' is replaced by '*'.
#' slot(aa, 'seqs')
#' 
#' 
#' aa <- peptide_phychem_index(base, acc = "MIYS850103", aaindex = "aaindex3")
#' aa
#' 
#' ## Description of the physicochemical index
#' slot(aa, 'phychem')
#' 
#' @author Robersy Sanchez <https://genomaths.com>
#' @export
setGeneric("peptide_phychem_index",
    function(
        aa,  
        ...) standardGeneric("peptide_phychem_index"))



## =========================== Characters ==========================

#' @aliases peptide_phychem_index
#' @rdname peptide_phychem_index
#' @param aa A character string, a \code{\link[Biostrings]{DNAStringSet}}
#' or a \code{\link[Biostrings]{DNAMultipleAlignment}} class object carrying
#' the DNA pairwise alignment of two sequences.
#' @param acc Accession id for a specified mutation or contact potential 
#' matrix.
#' @param aaindex Database where the requested accession id is locate and from
#' where the aminoacid indices can be obtained. The possible values are:
#' "aaindex2" or "aaindex3".
#' @param userindex User provided aminoacid indices. This can be a numerical
#' vector or a matrix (20 x 20). If a numerical matrix is provided, then the 
#' aminoacid indices are computes as column averages. 
#' @param alphabet Whether the alphabet is from the 20 aminoacid (AA) or
#' four (DNA)/RNA base alphabet. This would prevent mistakes, i.e., 
#' the strings "ACG" would be a base-triplet on the DNA alphabet or simply
#' the amino acid sequence of alanine, cysteine, and glutamic acid.
#' @param genetic.code,no.init.codon,if.fuzzy.codon The same as given in 
#' function [translation].
#' @return Depending on the user specifications, a mutation or contact 
#' potential matrix, a list of available matrices (indices) ids or index 
#' names can be returned. More specifically:
#' 
#' \describe{
#'  \item{\strong{aa_mutmat}: }{Returns an aminoacid mutation matrix or
#'    a statistical protein contact potentials matrix.}
#'  \item{\strong{aa_index}: }{Returns the specified aminoacid physicochemical 
#'    indices.}
#' }
#' 
#' @export
setMethod("peptide_phychem_index", signature(aa = "character"),
    function(
        aa, 
        acc = NULL,
        aaindex = NA,
        userindex = NULL,
        alphabet = c("AA", "DNA"),
        genetic.code = getGeneticCode("1"),
        no.init.codon = FALSE,
        if.fuzzy.codon = "error",
        ...) {
        
        alphabet <- match.arg(alphabet)
        
        AA <- c("A","R","N","D","C","Q","E","G","H",
                "I","L","K","M","F","P","S","T","W","Y","V")
    
        if (alphabet == "DNA") {
            nc <- all(nchar(aa) != 3)
            
            if (length(aa) == 1 && nchar(aa) > 1 &&  nc) 
                aa <- base2codon(aa)
            
            if (all(nchar(aa) != 3))
                stop("*** The argument 'a' is not a DNA protein-coding",
                    " sequence, i.e., it is not a base-triplet sequence.")
            
            aa <- translation(
                            x = aa,
                            genetic.code = genetic.code,
                            no.init.codon = no.init.codon,
                            if.fuzzy.codon = if.fuzzy.codon)
            
            alphabet <- "AA"
        }

        if (length(aa) == 1 && nchar(aa) > 1)
            aa <- str2chr(aa)
        
        if (!is.null(userindex)) {
            if (is.matrix(userindex)) 
                phychem <- colMeans(userindex)
            
            phychem <- phychem[ match(aa, AA) ]
        }
        else {
            if (is.null(acc))
                stop("The accesion ID for the aaindex must be provided")
            phychem <- aa_phychem_index(acc = acc, aaindex = aaindex)
            if (is.matrix(phychem))
                phychem <- colMeans(phychem)
            
            phychem <- phychem[ match(aa, AA) ]
        }
        return(phychem)
    }
)


## ========================== DNAStringSet ========================


#' @aliases peptide_phychem_index
#' @rdname peptide_phychem_index
#' @param num.cores,tasks Parameters for parallel computation using package
#' \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#' use, i.e. at most how many child processes will be run simultaneously
#' (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#' (only for Linux OS).
#' @param verbose If TRUE, prints the function log to stdout.
#' @param ... Not in use.
#' 
#' @export
setMethod("peptide_phychem_index", 
        signature(aa = "DNAStringSet_OR_DNAMultipleAlignment"),
    function(
        aa, 
        acc = NULL,
        aaindex = NA,
        userindex = NULL,
        alphabet = c("AA", "DNA"),
        genetic.code = getGeneticCode("1"),
        no.init.codon = FALSE,
        if.fuzzy.codon = "error",
        num.cores = 1L,
        tasks = 0L,
        verbose = FALSE,
        ...) {
        
        if (inherits(aa, "DNAMultipleAlignment"))
            aa <- aa@unmasked
        
        aa <- translation(
                        x = aa,
                        genetic.code = genetic.code,
                        no.init.codon = no.init.codon,
                        if.fuzzy.codon = if.fuzzy.codon)
        
        l <- width(aa[1])
        seqs <- as.character(aa)
        
        phychem <- character()
        if (is.na(aaindex))  
            aaindex <- "aaindex1"
        
        if (is.element(aaindex, c("aaindex1", "aaindex2", "aaindex3"))) {
            phychem <- aa_phychem_index(aaindex = aaindex, acc_list = TRUE)
            i <- grep(acc, phychem)
            if (length(i) > 0)
                phychem <- phychem[ i ]
            else
                stop("*** The accession provided '", acc, "' is not found in",
                    " the database '", aaindex, "'")
        }
        
        aa <- lapply(
                    seqs,
                    peptide_phychem_index, 
                    acc = acc, 
                    aaindex = aaindex,
                    userindex = userindex)
        
        
        nms <- names(seqs)
        aa <- matrix(unlist(aa, use.names = FALSE), 
                    ncol = l, byrow = TRUE)
        rownames(aa) <- paste0("S", seq(nrow(aa)))
        colnames(aa) <- paste0("A", seq(ncol(aa)))

        aa <- MatrixSeq(seqs = seqs,
                        matrix = aa, 
                        names = nms,
                        aaindex = if (is.na(aaindex)) "aaindex1" else aaindex,
                        phychem = phychem,
                        accession = acc)
        
        return(aa)
    }
)



