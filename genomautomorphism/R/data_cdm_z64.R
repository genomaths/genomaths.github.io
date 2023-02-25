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

#' Codon Distance Matrices for the Standard Genetic Code on Z4
#'
#' This is a list of 24 codon distance matrices created with function
#' \code{\link{codon_dist_matrix}} in the set of 24 genetic-code cubes on Z4
#' (using the default weights and assuming the standard genetic code (SGC). 
#' The data set is created to speed up the computation when working with DNA
#' sequences from superior organisms. Since distance matrices are symmetric, 
#' it is enough to provide the lower matrix. Each matrix is given as 
#' named/labeled vector (see the example).
#'
#' @examples 
#' ## Load the data set
#' data("cdm_z64", package = "GenomAutomorphism")
#' 
#' ## The lower matrix (given as vector) for cube "TCGA" (picking out the 20 
#' ## first values). Observe that this vector is labeled. Each numerical value 
#' ## corresponds to the distance between the codons specified by the 
#' ## name/label on it. For example, the distance between codons TTT and TCT 
#' ## is: 0.0625.
#' head(cdm_z64[[ "TCGA" ]], 20)
#' 
#' @format A list object.
"cdm_z64"
