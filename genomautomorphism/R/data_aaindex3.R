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

#' Statistical protein contact potentials matrices from AAindex ver.9.2
#' 
#' A statistical potential (also knowledge-based potential, empirical
#' potential, or residue contact potential) is an energy function derived from
#' an analysis of known structures in the Protein Data Bank.
#'
#' A list of 47 amino acid matrices from Amino Acid Index Database
#' \url{https://www.genome.jp/aaindex/} are provided here. AAindex is a
#' database of numerical indices representing various physicochemical and
#' biochemical properties of amino acids and pairs of amino acids.
#'
#' The contact potential matrix of amino acids is a set of 210 numerical
#' values, 20 diagonal and 20x19/2 off-diagonal elements, used for sequence
#' alignments and similarity searches.
#' 
#' @seealso [aaindex1], [aaindex2], and [get_mutscore].
#' @author Robersy Sanchez <https://genomaths.com>
#' @format A list carrying the the description 47 Amino Acid Matrices in 
#' AAindex ver.9.2 and the text file of matrices imported from
#' \url{https://www.genome.jp/aaindex/}.
#' @usage 
#' data("aaindex3", package = "GenomAutomorphism")
#' @examples 
#' ## Load the mutation matrices from database from the packages
#' data("aaindex3", package = "GenomAutomorphism")
#' 
#' ## Get the available mutation matrices
#' mat <- aa_mutmat(aaindex = "aaindex3", acc_list = TRUE)
#' mat[1:10]
#' 
"aaindex3"
