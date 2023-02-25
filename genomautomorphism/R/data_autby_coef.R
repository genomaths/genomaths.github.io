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
#' This is a \code{\link{AutomorphismList}} object carrying a list of pairwise
#' automorphisms between the DNA sequences from the MSA of primate somatic
#' cytochrome C grouped by automorphism's coefficients. The grouping derives
#' from the dataset \code{\link{brca1_autm}} after applying function
#' \code{\link{automorphism_bycoef}}.
#'
#' @format \code{\link{AutomorphismByCoefList}} class object.
"autby_coef"
