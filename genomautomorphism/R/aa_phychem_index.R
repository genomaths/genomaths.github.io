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

## =========================  Definition =======================

#' @rdname aa_phychem_index
#' @title Amino acid mutation matrix
#' @description 
#' The aminoacid similarity matrices from Amino Acid Index Database
#' \url{https://www.genome.jp/aaindex/} are provided here. 
#' \href{https://www.genome.jp/aaindex/}{AAindex} (ver.9.2) is a database of
#' numerical indices representing various physicochemical and biochemical
#' properties of amino acids and pairs of amino acids.
#'
#' The similarity of amino acids can be represented numerically, expressed in
#' terms of observed mutation rate or physicochemical properties. A similarity
#' matrix, also called a mutation matrix, is a set of 210 numerical values, 20
#' diagonal and 20x19/2 off-diagonal elements, used for sequence alignments 
#' and similarity searches.
#' 
#' Function \strong{\emph{aa_phychem_index}} is wrapper function to call two
#' other functions: \strong{\emph{aa_mutmat}} and \strong{\emph{aa_index}}
#' 
#' @param acc Accession id for a specified mutation or contact potential 
#' matrix.
#' @param aaindex Database where the requested accession id is locate. The 
#' possible values are:  "aaindex2" or "aaindex3".
#' @param acc_list Logical. If TRUE, then the list of available matrices ids 
#' and index names is returned.
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
#' @examples 
#' ## Load the mutation matrices from database from the packages
#' data("aaindex1","aaindex2", package = "GenomAutomorphism" )
#' 
#' ## Get the available mutation matrices
#' mat <- aa_mutmat(aaindex = "aaindex2", acc_list = TRUE)
#' mat[seq(10)]
#' 
#' ## Return the 'Base-substitution-protein-stability matrix 
#' ## (Miyazawa-Jernigan, 1993)'
#' aa_mutmat(acc = "MIYS930101", aaindex = "aaindex2")
#' 
#' ## Return the 'BLOSUM80 substitution matrix (Henikoff-Henikoff, 1992)'
#' aa_mutmat(acc = "HENS920103", aaindex = "aaindex2")
#' 
#' ## Using wrapping function
#' aa_phychem_index(acc = "EISD840101", aaindex = "aaindex1")
#' 
#' ## Just the info. The information provided after the reference
#' ## corresponds to the correlaiton of 'EISD840101' with other indices.
#' aa_phychem_index(acc = "EISD840101", aaindex = "aaindex1", info = TRUE)
#' 
#' @seealso [aaindex1], [aaindex2],  [aaindex3], and [get_mutscore].
#' @author Robersy Sanchez <https://genomaths.com>
#' @aliases aa_phychem_index
#' @export
#' 
aa_phychem_index <- function(
        acc = NA, 
        aaindex = NA,
        acc_list = FALSE,
        info = FALSE) {
    
    if (is.na(aaindex) || aaindex == "aaindex1")
        res <- aa_index(
            acc = acc, 
            acc_list = acc_list,
            info = info)
    else {
        res <- aa_mutmat(
            acc = acc, 
            aaindex = aaindex,
            acc_list = acc_list)
    }
    
    return(res)
}
    

## ===================== aa_mutmat ====================


#' @rdname aa_phychem_index
#' @aliases aa_mutmat
#' @export
aa_mutmat <- function(
    acc = NA, 
    aaindex = c("aaindex2", "aaindex3"),
    acc_list = FALSE) {
    
    aaindex <- match.arg(aaindex)
    switch(
        aaindex,
        "aaindex2" = data("aaindex2", package = "GenomAutomorphism",
                        envir = environment()),
        "aaindex3" = data("aaindex3", package = "GenomAutomorphism",
                        envir = environment())
    )
    
    aaindex <- get(aaindex)
    
    if (acc_list)
        return(aaindex$acc_num)
    else {
        if (is.na(acc))
            stop("*** Please provide an accession id: 'acc'.",
                " See the function examples.")
        
        idx <- grep(acc, aaindex$aaindex)
        
        if (length(idx) == 0)
            stop("*** The accession provided '", acc, "' is not found in",
                " the database '", aaindex, "'")
            
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

## ===================== aa_index ====================

#' @rdname aa_phychem_index
#' @aliases aa_index
#' @param info Logical. if TRUE, then whole information for the 
#' physicochemical index will be returned.
#' @export
aa_index <- function(
        acc = NA, 
        acc_list = FALSE,
        info = FALSE) {
    
    data("aaindex1", package = "GenomAutomorphism",
        envir = environment())
    aaindex1 <- get("aaindex1")
    
    if (acc_list)
        return(aaindex1$acc_num)
    else {
        if (is.na(acc))
            stop("*** Please provide an accession id: 'acc'.",
                " See the function examples.")
        
        i0 <- idx <- match(paste0("H ", acc), aaindex1$aaindex)
        if (is.na(i0))
            stop("*** The accession provided '", acc, "' is not found in",
                " the database 'aaindex1'.")
            
        patt <- FALSE
        while (!patt) {
            idx <- idx + 1
            patt <- grepl("^I", aaindex1$aaindex[idx])
        }
        
        if (info)
            return(aaindex1$aaindex[seq(i0, idx + 2, 1)])
        
        nums <- aaindex1$aaindex[seq(idx + 1, idx + 2, 1)]
        nums <- gsub("   ", ",", nums)
        nums <- slapply(nums, function(s) {
            s <- na.omit(as.numeric(strsplit(s, ",")[[1]]))
        })
        names(nums) <- c("A","R","N","D","C","Q","E","G","H",
                        "I","L","K","M","F","P","S","T","W","Y","V")
        return(nums)
    }
}

