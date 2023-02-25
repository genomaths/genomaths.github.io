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

#' List of 94 Amino Acid Matrices from AAindex
#'
#' The aminoacid similarity matrices from Amino Acid Index Database
#' \url{https://www.genome.jp/aaindex/} are provided here. AAindex (ver.9.2) 
#' is a database of numerical indices representing various physicochemical and
#' biochemical properties of amino acids and pairs of amino acids.
#'
#' The similarity of amino acids can be represented numerically, expressed in
#' terms of observed mutation rate or physicochemical properties. A similarity
#' matrix, also called a mutation matrix, is a set of 210 numerical values, 20
#' diagonal and 20x19/2 off-diagonal elements, used for sequence alignments 
#' and similarity searches.

#' @examples 
#' ## Load the mutation matrices from database from the packages
#' data(aaindex2, package = "GenomAutomorphism")
#' 
#' ## Get the available mutation matrices
#' mat <- aa_mutmat(aaindex = aaindex2, acc_list = TRUE)
#' mat[1:10]
#' @seealso \code{\link{aaindex2}} and \code{\link{aa_mutmat}}, and
#' \code{\link{get_mutscore}}.
#' @author Robersy Sanchez <https://genomaths.com>
#' @format \code{\link{AutomorphismList}} class object.
"aaindex2"
