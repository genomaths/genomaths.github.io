## Copyright (C) 2021-2022 Robersy Sanchez <https://genomaths.com/>
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

#' Replace bases with integers from Z4 and Z5
#' @rdname base2int
#' @aliases base2int
#' @description A simple function to represent DNA bases as elements from 
#' the Abelian group of integers modulo 4 (Z4) or 5 (Z5).
#' @param base A character vector, string , or a dataframe of letters from the 
#' DNA/RNA alphabet.
#' @param group A character string denoting the group representation for the
#' given base or codon as shown in reference (2-3).
#' @param cube A character string denoting one of the 24 Genetic-code cubes,
#' as given in references (2-3).
#' @param ... Not in use.
#' @return A numerical vector.
#' @aliases base2int
#' @author Robersy Sanchez <https://genomaths.com>
#' @seealso \code{\link{base_coord}} and \code{\link{codon_coord}}.
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
#' ## data.frames must carry only single letters
#' \donttest{
#' base2int(data.frame(x1 = c("UTG", "GTA"), x2 = c("UTG", "GTA")))
#' }
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

#' @rdname base2int
#' @aliases base2int
#' @export
setMethod("base2int", signature(base = "character"),
    function(
        base, 
        group = c("Z4", "Z5", "Z64", "Z125", "Z4^3", "Z5^3"),
        cube = c(
            "ACGT", "AGCT", "TCGA", "TGCA", "CATG",
            "GTAC", "CTAG", "GATC", "ACTG", "ATCG",
            "GTCA", "GCTA", "CAGT", "TAGC", "TGAC",
            "CGAT", "AGTC", "ATGC", "CGTA", "CTGA",
            "GACT", "GCAT", "TACG", "TCAG")) {
        
        cube <- match.arg(cube)
        group <- match.arg(group)
        
        if (length(base) > 1) { 
            base <- str2chr(base)
            base <- lapply(base, base_repl, cube = cube, group = group)
            base <- do.call(rbind, base)
        }
        else {
            base <- str2chr(base)
            base <- base_repl(base = base, cube = cube, group = group)
        }
        return(base)
    }
)


#' @rdname base2int
#' @aliases base2int
#' @export
setMethod("base2int", signature(base = "data.frame"),
    function(
        base,
        group = c("Z4", "Z5", "Z64", "Z125", "Z4^3", "Z5^3"), 
        cube = c(
            "ACGT", "AGCT", "TCGA", "TGCA", "CATG",
            "GTAC", "CTAG", "GATC", "ACTG", "ATCG",
            "GTCA", "GCTA", "CAGT", "TAGC", "TGAC",
            "CGAT", "AGTC", "ATGC", "CGTA", "CTGA",
            "GACT", "GCAT", "TACG", "TCAG")) {
        
        if (any(nchar(as.matrix(base)) > 1))
            stop("*** Argument 'x' must carry only single letters.")
        
        cube <- match.arg(cube)
        group <- match.arg(group)
        base <- base_repl(base = base, cube = cube, group = group)
        return(base)
    }
)



## --------------------------- Auxiliary functions --------------------------

#' Replace bases with integers
#' @details Internal use only.
#' @keywords internal
#' @return A numerical vector.
base_repl <- function(base, cube, group) {
    alf <- strsplit(cube, "")[[1]]
    
    if (group == "Z4") {
        base[base == "U"] <- "T"
        base[base == alf[1]] <- 0
        base[base == alf[2]] <- 1
        base[base == alf[3]] <- 2
        base[base == alf[4]] <- 3
        base[base == "-"] <- NA
        base[base == "N"] <- NA
    }
    
    if (group == "Z5") {
        base[base == "U"] <- "T"
        base[base == alf[1]] <- 1
        base[base == alf[2]] <- 2
        base[base == alf[3]] <- 3
        base[base == alf[4]] <- 4
        base[base == "-"] <- 0
        base[base == "N"] <- 0
        base[base == "D"] <- 0
    }
    
    if (is.data.frame(base)) {
        base <- apply(base, 2, as.numeric)
    } else {
        base <- as.numeric(base)
    }
    return(base)
}


