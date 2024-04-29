## Copyright (C) 2022-2024 Robersy Sanchez <https://genomaths.com/>
## Author: Robersy Sanchez This file is part of the R package
## 'GenomAutomorphism'.  'GenomAutomorphism' is a free
## software: you can redistribute it and/or modify it under the
## terms of the GNU General Public License as published by the Free
## Software Foundation, either version 3 of the License, or (at
## your option) any later version.  This program is distributed in
## the hope that it will be useful, but WITHOUT ANY WARRANTY;
## without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR P


## =========================== Definition ==========================

#' @rdname dna_phychem
#' @aliases dna_phychem
#' @title DNA numerical matrix
#' @param seqs A character string, a \code{\link[Biostrings]{DNAStringSet}}
#' or a \code{\link[Biostrings]{DNAMultipleAlignment}} class object carrying
#' the DNA pairwise alignment of two sequences.
#' @param phychem A list of DNA bases physicochemical properties, e.g., like
#' those provided in [dna_phyche].
#' @param index_name Optional. Name of breve description of the base
#' physicochemical property applied to represent the DNA sequence.
#' @param ... Not in use.
#' @description 
#' This function applies the numerical indices representing various
#' physicochemical and biochemical properties of DNA bases. As results, DNA
#' sequences are represented as numerical vectors which can be subject of
#' further downstream statistical analysis and digital signal processing.
#' @author Robersy Sanchez <https://genomaths.com>
#' @seealso [peptide_phychem_index]
#' @returns A [MatrixSeq]-class object.
#' @export
#' @examples
#' ## Let's create DNAStringSet-class object 
#' base <- DNAStringSet(x = c( seq1 ='ACGTGATCAAGT', 
#'                             seq2 = 'GTGTGATCCAGT', 
#'                             seq3 = 'TCCTGATCAGGT'))
#' 
#' 
#' dna_phychem(seqs = base,
#'             phychem = list('A' = 0.87, 'C' = 0.88, 'T' = 0.82,
#'                             'G' = 0.89, 'N' = NA),
#'             index_name = "Proton-Affinity")
#' 
setGeneric("dna_phychem",
    function(
        seqs,  
        ...) standardGeneric("dna_phychem"))

## ====================== Character Method  =======================

#' @aliases dna_phychem
#' @rdname dna_phychem
#' @export
setMethod("dna_phychem", signature(seqs = "character"), 
    function(
        seqs,
        phychem = list('A' = NULL, 'T' = NULL, 'C' = NULL,
                        'G' = NULL, 'N' = NULL)) {
        
        base2int(base = seqs, phychem = phychem)
    }
)


## ========================== DNAStringSet ========================


#' @aliases dna_phychem
#' @rdname dna_phychem
#' @export
setMethod("dna_phychem", 
    signature(seqs = "DNAStringSet_OR_DNAMultipleAlignment"),
    function(
        seqs, 
        phychem = list('A' = NULL, 'T' = NULL, 'C' = NULL,
                        'G' = NULL, 'N' = NULL),
        index_name = NULL,
        ...) {
        
        if (inherits(seqs, "DNAMultipleAlignment"))
            seqs <- seqs@unmasked
        
        nms <- names(seqs)
        
        seqs <- as.character(seqs)
        base <- base2int(base = seqs, phychem = phychem)
        
        rownames(base) <- paste0("S", seq(nrow(base)))
        colnames(base) <- paste0("B", seq(ncol(base)))
        
        if (is.null(index_name))
            index_name <- character()
        
        base <- MatrixSeq(
                        seqs = seqs,
                        matrix = base, 
                        names = nms,
                        aaindex = character(),
                        phychem = index_name,
                        accession = character())
        
        return(base)
    }
)


