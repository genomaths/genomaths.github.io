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

#' @rdname translation
#' @title Translation of DNA/RNA sequences
#' @description This function extends \code{\link[Biostrings]{translate}} 
#' function to include letters that are frequently found in the DNA sequence
#' databases to indicate missing information and are not part of the the 
#' DNA/RNA alphabet. Also, it is able to process sequences as just simple 
#' 'character' objects. 
#' @details If argument 'x' belong to any of the classes admitted by function 
#' \code{\link[Biostrings]{translate}}, then this function is called to make
#' the translation.
#' @param x A character string or the same arguments given to function
#' \code{\link[Biostrings]{translate}}.
#' @param genetic.code The same as in \code{\link[Biostrings]{translate}}
#' @param no.init.codon,if.fuzzy.codon Used only if 'x' is not a 'character'
#' object. The same as in \code{\link[Biostrings]{translate}}.
#' @param ... Not in use yet.
#' @seealso \code{\link[Biostrings]{translate}}
#' @author Robersy Sanchez <https://genomaths.com>
#' @import Biostrings
#' @examples 
#' ## Load a small DNA sequence alingment 
#' data("aln", package = "GenomAutomorphism")
#' 
#' translation(aln)
#' 
#' ## Load a pairwise DNA sequence alingment of COVID-19 genomes
#' data("covid_aln", package = "GenomAutomorphism")
#' 
#' translation(covid_aln)
#' @return The translated amino acid sequence.
#' @export
setGeneric(
    "translation",
    function(x,
            ...) {
        standardGeneric("translation")
    }
)

#' @rdname translation
#' @aliases translation
#' @export
setMethod(
    "translation", signature(x = "character"),
    function(
        x,
        genetic.code = getGeneticCode("1"),
        no.init.codon = FALSE,
        if.fuzzy.codon = "error") {
        
        if (length(x) > 1 && any((nchar(x) > 1))) {
            x <- slapply(x, function(s) paste0(transl(s), collapse = ""))
        }else
            x <- transl(x)
        return(x)
    }
)

setClassUnion(
    "BioString", 
    c("DNAStringSet","RNAStringSet","DNAStringSet", 
        "RNAStringSet", "DNAString", "RNAString", "MaskedDNAString",
        "MaskedRNAString", "DNAString", "RNAString", "MaskedDNAString", 
        "MaskedRNAString","DNAMultipleAlignment"))


#' @rdname translation
#' @aliases translation
#' @export
setMethod(
    "translation", signature(x = "BioString"),
    function(
        x,
        genetic.code = getGeneticCode("1"),
        no.init.codon = FALSE,
        if.fuzzy.codon = "error") {
        
        if (inherits(x, "DNAMultipleAlignment")) 
            x <- unmasked(x)
        
        aa <- try(translate(
            x,
            genetic.code = genetic.code,
            no.init.codon = no.init.codon,
            if.fuzzy.codon = if.fuzzy.codon), 
            silent = TRUE)
        
        if (inherits(aa, "try-error")) {
            aa <- translation(as.character(x))
            rm(x); gc()
            aa <- AAStringSet(aa)
        }
        return(aa)
    }
)

## ========== Auxiliary function ============

transl <- function(
        x, 
        genetic.code = getGeneticCode("1")) {
    
    if (length(x) == 1) {
        if (nchar(x) %% 3 != 0)
            stop(
                "*** Argument 'x' must be a character vector or ",
                "a character coercible to a character vector"
            )
        x <- base2codon(x)
    }
    
    if (length(x) > 1) {
        if (all(nchar(x) == 1)) {
            x <- paste(x, collapse = "")
            x <- base2codon(x)
        }
    }
    
    if (any((nchar(x) %% 3) != 0))
        stop("*** The number of characters in argument 'x' ,
                    must be multiple of 3")
    
    x <- toupper(x)
    x <- gsub("U", "T", x)
    
    x <- genetic.code[match(x, names(genetic.code))]
    x[is.na(x)] <- "-"
    return(x)
}

