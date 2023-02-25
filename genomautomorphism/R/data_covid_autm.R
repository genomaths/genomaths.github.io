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

#' Automorphisms between DNA Sequences from two COVID-19 genomes
#'
#' This is a \code{\link{AutomorphismList}} object carrying a list of pairwise
#' automorphisms between the SARS coronavirus GZ02 (GenBank: AY390556.1:
#' 265-13398_13398-21485) and Bat SARS-like coronavirus isolate bat-SL-CoVZC45
#' (GenBank: KY417151.1: protein-coding regions). The pairwise
#' DNA sequence alignment is available in the dataset named
#' \code{\link{covid_aln}} and the automorphisms were estimated with function
#' \code{\link{autZ64}}.
#'
#' @format \code{\link{AutomorphismList}} class object.
"covid_autm"
