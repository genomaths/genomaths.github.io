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
#' This is a \code{\link{AutomorphismList}} object carrying a list of pairwise
#' automorphisms between the DNA sequences from the MSA of primate BRCA1
#' DNA repair gene. The data set brca1_aln2 has 41 DNA sequences and it 
#' contains the previous 20 primate variants found in 'braca1_aln' data set
#' plus 21 single mutation variants (SMV) from the human sequence NM_007298 
#' transcript variant 4. The location of each SMV is given in the heading from
#' each sequence.
#' 
#' The automorphisms were estimated from the \code{\link{brca1_aln}} MSA with
#' function \code{\link{autZ64}}.
#'
#' @format \code{\link{AutomorphismList}} class object.
"brca1_autm2"
