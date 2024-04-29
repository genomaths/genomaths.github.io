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

#' List of 571 Amino Acid Physicochemical Indexes from AAindex Database
#'
#' The aminoacid indexes from Amino Acid Index Database
#' \url{https://www.genome.jp/aaindex/} are provided here. AAindex (ver.9.2) 
#' is a database of numerical indices representing various physicochemical and
#' biochemical properties of amino acids and pairs of amino acids.
#' 
#' @seealso [aaindex2] and [aaindex3].
#' @author Robersy Sanchez <https://genomaths.com>
#' @format A list carrying the the description 566 Amino Acid Indices in 
#' AAindex ver.9.2 and the text file with the matrices imported from
#' \url{https://www.genome.jp/aaindex/}.
#' @docType data
#' @usage 
#' data("aaindex1", package = "GenomAutomorphism")
#' @examples 
#' ## Load the mutation matrices from database from the packages
#' data("aaindex1", package = "GenomAutomorphism", envir = environment())
#' 
#' ## Get the available aminoacid indices.
#' mat <- aa_phychem_index(aaindex = "aaindex1", acc_list = TRUE)
#' mat[1:10]
#' 
"aaindex1"
