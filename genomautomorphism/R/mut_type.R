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

#' Classification of DNA base mutations
#' @description Each DNA/RNA base can be classified into three main classes
#' according to three criteria (1): number of hydrogen bonds (strong-weak),
#' chemical type (purine-pyrimidine), and chemical groups (amino versus keto).
#' Each criterion produces a partition of the set of bases: 1) According to the
#' number of hydrogen bonds (on DNA/RNA double helix): strong S={C,G} (three
#' hydrogen bonds) and weak W={A,U} (two hydrogen bonds). According to the
#' chemical type: purines R={A, G} and pyrimidines Y={C,U}. 3). According to
#' the presence of amino or keto groups on the base rings: amino M={C,A} and
#' keto K={G,U}. So, each mutational event can be classified as according to
#' the type of involved in it (2).

#' @param x,y Character strings denoting DNA bases
#' @return A character string of same length of 'x' and 'y'.
#' @export
#' @examples
#' ## Mutation type 'R'
#' mut_type("A", "G")
#' 
#' ## Mutation type 'M'
#' mut_type("A", "C")
#' 
#' ## Mutation type 'W'
#' mut_type("A", "T")
#' 
#' ## Mutation type 'S'
#' mut_type("G", "C")
#' @references
#' \enumerate{
#'      \item A. Cornish-Bowden, Nomenclature for incompletely specified bases
#'          in nucleic acid sequences: recommendations 1984, Nucleic Acids
#'          Res. 13 (1985) 3021-3030.
#'      \item MA.A. Jimenez-Montano, C.R. de la Mora-Basanez, T. Poschel, The
#'          hypercube structure of the genetic code explains conservative and
#'          non-conservative aminoacid substitutions in vivo and in vitro,
#'          Biosystems. 39 (1996) 117-125.
#' }
#'
mut_type <- function(x, y) {
    if (nchar(x) != nchar(y)) {
        stop("*** Arguments 'x' & 'y' must have the same length.")
    }
    res <- mapply(function(b1, b2) {
        if (b1 != b2) {
            bp <- paste0(b1, b2)
            idx <- match(bp, mut_type1)
            if (is.na(idx)) {
                idx <- match(bp, mut_type2)
                res <- names(mut_type2[idx])
            } else {
                res <- names(mut_type1[idx])
            }
        } else {
            res <- "H"
        }
        return(res)
    }, str2chr(x), str2chr(y))
    res <- paste(res, collapse = "")
    if (grepl("NA", res)) {
        res <- gsub("NA", "-", res)
    }
    return(res)
}


## =========== Auxiliary function & data ======================

# str2ch <- function(x) {
#     strsplit(x, split = "")[[1]]
# }

mut_type1 <- c(
    R = paste0("A", "G"),
    Y = paste0("C", "T"),
    S = paste0("C", "G"),
    W = paste0("A", "T"),
    M = paste0("A", "C"),
    K = paste0("G", "T")
)

mut_type2 <- c(
    R = paste0("G", "A"),
    Y = paste0("T", "C"),
    S = paste0("G", "C"),
    W = paste0("T", "A"),
    M = paste0("C", "A"),
    K = paste0("T", "G")
)
