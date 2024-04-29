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

#' Automorphisms between DNA Sequences from Primate BRCA1 Genes
#'
#' This is a [AutomorphismList] object carrying a list of pairwise
#' automorphisms between the DNA sequences from the MSA of primate BRCA1
#' DNA repair gene. The automorphisms were estimated from the
#' [brca1_aln] MSA with function [autZ64].
#' @usage 
#' data("brca1_autm", package = "GenomAutomorphism")
#'
#' @format [AutomorphismList] class object.
#' @author Robersy Sanchez <https://genomaths.com>
#' @seealso [brca1_autm2], [brca1_aln], [brca1_aln2], and [covid_autm].
#' @examples
#' data("brca1_autm", package = "GenomAutomorphism")
#' brca1_autm
#' 
"brca1_autm"
