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

## ===================== Definition ==========================


#' Replace bases with integers from Z4 and Z5
#' @rdname base2int
#' @aliases base2int
#' @description A simple function to represent DNA bases as elements from 
#' the Abelian group of integers modulo 4 (Z4), 5 (Z5), or 2 (Z2).
#' @details
#' For Z2 (binary representation of DNA bases), the cube bases are represented
#' in their order by: '00', '01', '10', and '11' (examples section). 
#' 
#' @param base A character vector, string , or a dataframe of letters from the 
#' DNA/RNA alphabet.
#' @param group A character string denoting the group representation for the
#' given base or codon as shown in reference (2-3).
#' @param cube A character string denoting one of the 24 Genetic-code cubes,
#' as given in references (2-3).
#' @param phychem Optional. Eventually, it could be useful to represent 
#' DNA bases by numerical values of measured physicochemical properties. If
#' provided, then this argument must be a named numerical list. For example,
#' the \code{\link[base]{scale}} values of deoxyribonucleic acids proton 
#' affinity (available at \url{https://www.wolframalpha.com/} and
#' in cell phone app: Wolfram Alpha):
#' 
#' \eqn{list('A' = 0.87, 'C' = 0.88, 'T' = 0.82, 'G' = 0.89, 'N' = NA)}
#' 
#' where symbol 'N' provide the value for any letter out of DNA base alphabet.
#' In this example, we could write NA or 0 (see example section).
#' 
#' @param ... Not in use.
#' @return A numerical vector.
#' @aliases base2int
#' @author Robersy Sanchez <https://genomaths.com>
#' @seealso [base_coord], [codon_coord], and [dna_phychem].
#' @examples 
#' ## A triplet with a letter not from DNA/RNA alphabet
#' ## 'NA' is introduced by coercion!
#' base2int("UDG")
#' 
#' ## The base replacement in cube "ACGT and group "Z4"
#' base2int("ACGT")
#' 
#' ## The base replacement in cube "ACGT and group "Z5"
#' base2int("ACGT", group = "Z5")
#' 
#' ## A vector of DNA base triplets
#' base2int(c("UTG", "GTA"))
#' 
#' ##  A vector of DNA base triplets with different number of triplets.
#' ##  Codon 'GTA' is recycled!
#' base2int(base = c("UTGGTA", "CGA"), group = "Z5")
#' 
#' ## Data frames 
#' 
#' base2int(data.frame(x1 = c("UTG", "GTA"), x2 = c("UTG", "GTA")))
#' 
#' 
#' ## Cube bases are represented n their order by: '00', '01', '10', and '11',
#' ## For example for cube = "ACGT" we have mapping: A -> '00', C -> '01',
#' ## G -> '11', and C -> '10'.
#' 
#' base2int("ACGT", group = "Z2", cube = "ACGT")
#' 
#' @references
#' \enumerate{
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
#' @export
setGeneric(
    "base2int",
    function(
        base,
        ...) {
        standardGeneric("base2int")
    }
)

## ===================== Character method ======================

#' @rdname base2int
#' @aliases base2int
#' @export
setMethod("base2int", signature(base = "character"),
    function(
        base, 
        group = c("Z4", "Z5", "Z64", "Z125", "Z4^3", "Z5^3", "Z2"),
        cube = c(
            "ACGT", "AGCT", "TCGA", "TGCA", "CATG",
            "GTAC", "CTAG", "GATC", "ACTG", "ATCG",
            "GTCA", "GCTA", "CAGT", "TAGC", "TGAC",
            "CGAT", "AGTC", "ATGC", "CGTA", "CTGA",
            "GACT", "GCAT", "TACG", "TCAG"),
        phychem = list('A' = NULL, 'T' = NULL, 'C' = NULL,
                    'G' = NULL, 'N' = NULL)) {
        
        cube <- match.arg(cube)
        group <- match.arg(group)
        
        if (length(base) > 1) { 
            base <- str2chr(base)
            base <- lapply(
                        base, 
                        base_repl, 
                        cube = cube, 
                        group = group,
                        phychem = phychem)
            base <- do.call(rbind, base)
        }
        else {
            base <- str2chr(base)
            base <- base_repl(
                            base = base, 
                            cube = cube, 
                            group = group,
                            phychem = phychem)
        }
        return(base)
    }
)

## ===================== Data frame method ======================

#' @rdname base2int
#' @aliases base2int
#' @export
setMethod("base2int", signature(base = "data.frame"),
    function(
        base,
        group = c("Z4", "Z5", "Z64", "Z125", "Z4^3", "Z5^3", "Z2"), 
        cube = c(
            "ACGT", "AGCT", "TCGA", "TGCA", "CATG",
            "GTAC", "CTAG", "GATC", "ACTG", "ATCG",
            "GTCA", "GCTA", "CAGT", "TAGC", "TGAC",
            "CGAT", "AGTC", "ATGC", "CGTA", "CTGA",
            "GACT", "GCAT", "TACG", "TCAG"),
        phychem = list('A' = NULL, 'T' = NULL,  'C' = NULL, 
                    'G' = NULL, 'N' = NULL)) {
        
        nch <- nchar(as.matrix(base))
        if (any(nchar(as.matrix(base)) > 1)) {
            cn <- colnames(base)
            rn <- rownames(base)
            
            nc <- ncol(base)
            nr <- nrow(base)
            
            nch <- unique(as.vector(nch))
            if (length(nch) > 1)
                stop("*** All the strings in argument 'base' must have",
                    "the same number of letters.")

            base <- try(apply(base, 1, str2chr ), silent = TRUE)
            if (inherits(base, "try-error"))
                stop("*** Argument 'base' must carry only single letters.",
                    "It was not possible to transform string into letters.")
            
            cn <- slapply(cn, function(x) rep(x, nch))
            base <- data.frame(
                        matrix(
                            unlist(base, use.names = FALSE), 
                                ncol = nch * nc))
            
            colnames(base) <- cn
            rownames(base) <- rn
            
        }
        
        cube <- match.arg(cube)
        group <- match.arg(group)
        base <- base_repl(
                        base = base, 
                        cube = cube, 
                        group = group,
                        phychem = phychem)
        return(base)
    }
)



## --------------------------- Auxiliary functions --------------------------

#' Replace bases with integers
#' @details Internal use only.
#' @keywords internal
#' @return A numerical vector.
base_repl <- function(
    base, 
    cube, 
    group,
    phychem) {
    
    z4 <- setdiff(toupper(letters), c("A", "C", "G", "T", "U"))
    
    if (all(slapply(phychem, is.null))) {
        
        alf <- strsplit(cube, "")[[1]]
        
        if (group == "Z4") {
            base[base == "U"] <- "T"
            base[base == alf[1]] <- 0
            base[base == alf[2]] <- 1
            base[base == alf[3]] <- 2
            base[base == alf[4]] <- 3
            base[base == "-"] <- NA
            base[base == "N"] <- NA
            base[base == "?"] <- NA
            base[base == "*"] <- NA
            base[is.element(base, z4)] <- NA
        }
        
        if (group == "Z5") {
            z5 <- setdiff(toupper(letters), c("A", "C", "D", "G", "T", "U"))
            base[base == "U"] <- "T"
            base[base == alf[1]] <- 1
            base[base == alf[2]] <- 2
            base[base == alf[3]] <- 3
            base[base == alf[4]] <- 4
            base[base == "-"] <- 0
            base[base == "N"] <- 0
            base[base == "D"] <- 0
            base[base == "?"] <- 0
            base[base == "*"] <- 0
            base[is.element(base, z5)] <- 0
        }
        
        if (group == "Z2") {
            base[base == "U"] <- "T"
            base[base == alf[1]] <- '00'
            base[base == alf[2]] <- '01'
            base[base == alf[3]] <- '10'
            base[base == alf[4]] <- '11'
            base[base == "-"] <- '00'
            base[base == "N"] <- '00'
            base[base == "D"] <- '00'
            base[base == "?"] <- '00'
            base[base == "*"] <- '00'
            base[is.element(base, z4)] <- '00'
            base <- as.integer(strsplit(paste0(base, collapse = ""),
                            split = "")[[1]])
        }
    }
    else {
        base[base == "U"] <- "T"
        base[base == "A"] <- phychem$A
        base[base == "T"] <- phychem$T
        base[base == "C"] <- phychem$C
        base[base == "G"] <- phychem$G
        base[base == "-"] <- phychem$N
        base[base == "N"] <- phychem$N
        base[base == "D"] <- phychem$N
        base[base == "?"] <- phychem$N
        base[base == "*"] <- phychem$N
        base[is.element(base, z4)] <- phychem$N
    }
    
    if (is.data.frame(base)) {
        base <- apply(base, 2, as.numeric)
    } else {
        base <- as.numeric(base)
    }
    return(base)
}


