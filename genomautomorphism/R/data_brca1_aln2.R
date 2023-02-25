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

#' Multiple Sequence Alignment (MSA) of Primate BRCA1 DNA repair genes.
#'
#' This is a \code{\link[Biostrings]{DNAMultipleAlignment}} carrying a MSA of
#' [BRCA1 DNA repair genes](https://bit.ly/3DimROD) to be used in the
#' examples provided for the package functions. The original file can be
#' downloaded from GitHub at: <https://bit.ly/3DimROD>. This data set has
#' 41 DNA sequences and it contains the previous 20 primate variants found in
#' 'brca1_aln' data set plus 21 single mutation variants (SMV) from the human
#' sequence NM_007298 transcript variant 4. The location of each SMV is given
#' in the heading from each sequence.
#'
#' @format \code{\link[Biostrings]{DNAMultipleAlignment}} class object.
"brca1_aln2"
