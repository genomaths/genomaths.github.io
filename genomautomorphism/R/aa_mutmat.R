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

#' @rdname aa_mutmat
#' @title Amino acid mutation matrix
#' @description This returns an amino acid mutation matrix or a
#' statistical protein contact potentials matrix from 
#' \href{https://www.genome.jp/aaindex/}{AAindex} (ver.9.2).
#' 
#' The aminoacid similarity matrices from Amino Acid Index Database
#' \url{https://www.genome.jp/aaindex/} are provided here. AAindex (ver.9.2) 
#' is a database of numerical indices representing various physicochemical and
#' biochemical properties of amino acids and pairs of amino acids.
#'
#' The similarity of amino acids can be represented numerically, expressed in
#' terms of observed mutation rate or physicochemical properties. A similarity
#' matrix, also called a mutation matrix, is a set of 210 numerical values, 20
#' diagonal and 20x19/2 off-diagonal elements, used for sequence alignments 
#' and similarity searches.
#' @param acc Accession id for a specified mutation or contact potential 
#' matrix.
#' @param aaindex Database where the requested accession id is locate. The 
#' possible values are:  "aaindex2" or "aaindex3".
#' @param acc_list Logical. If TRUE, then the list of available matrices ids 
#' and index names is returned.
#' @return A mutation or contact potential matrix, or the list of available
#' matrices ids and index names is returned.
#' @examples 
#' ## Load the mutation matrices from database from the packages
#' data("aaindex2", package = "GenomAutomorphism" )
#' 
#' ## Get the available mutation matrices
#' mat <- aa_mutmat(aaindex = aaindex2, acc_list = TRUE)
#' mat[1:10]
#' 
#' ## Return the 'Base-substitution-protein-stability matrix 
#' ## (Miyazawa-Jernigan, 1993)'
#' aa_mutmat(acc = "MIYS930101", aaindex = aaindex2)
#' 
#' ## Return the 'BLOSUM80 substitution matrix (Henikoff-Henikoff, 1992)'
#' aa_mutmat(acc = "HENS920103", aaindex = aaindex2)
#' 
#' @seealso \code{\link{aaindex2}}, \code{\link{aaindex3}}, and
#' \code{\link{get_mutscore}}.
#' @author Robersy Sanchez <https://genomaths.com>
#' @export
aa_mutmat <- function(
    acc = NA, 
    aaindex = NA,
    acc_list = FALSE) {
    
    if (acc_list)
        return(aaindex$acc_num)
    else {
        idx <- grep(acc, aaindex$aaindex)
        patt <- FALSE
        while (!patt) {
            idx <- idx + 1
            patt <- grepl("^M", aaindex$aaindex[idx])
        }
        
        nums <- aaindex$aaindex[seq(idx + 1, idx + 20, 1)]
        
        nums <- lapply(seq(nums), function(k) {
            s <- strsplit(nums[k], " ")[[1]]
            s <- as.numeric(s[nchar(s) > 0])
            s <- c(s, rep(NA, 20 - length(s)))
            return(s)
        })
        rown <- c("A","R","N","D","C","Q","E","G","H",
                "I","L","K","M","F","P","S","T","W","Y","V")
        
        nums <- do.call(rbind, nums)
        rownames(nums) <- rown
        colnames(nums) <- rown
        nums[upper.tri(nums)] <- t(nums)[upper.tri(nums)]
        return(nums)
    }
}
