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

#' Automorphisms between DNA Primate BRCA1 Genes Grouped by Coefficients
#'
#' This is a [AutomorphismList] object carrying a list of pairwise
#' automorphisms between the DNA sequences from the MSA of primate somatic
#' cytochrome C grouped by automorphism's coefficients. The grouping derives
#' from the dataset [brca1_autm] after applying function
#' [automorphism_bycoef].
#'
#' @format [AutomorphismByCoefList] class object.
#' @usage 
#' data("autby_coef", package = "GenomAutomorphism")
#' @examples
#' ## Load the data set
#' data("autby_coef", package = "GenomAutomorphism")
#' autby_coef
#' 
#' ## Mutation type found in the data
#' unique(autby_coef$human_1.human_2$mut_type)
#' 
"autby_coef"
