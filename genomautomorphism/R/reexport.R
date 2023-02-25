#' Reexport useful functions to be available to users
## From S4Vectors ---------------------------------------

#' @importFrom S4Vectors mcols
#' @export
#' @examples 
#' ## Load an Automorphism object and take its metacolumns
#' data("autm", package = "GenomAutomorphism")
#' mcols(autm)
S4Vectors::mcols

#' @importFrom S4Vectors mcols<-
#' @export
S4Vectors::`mcols<-`

#' @importFrom S4Vectors setValidity2
#' @export
S4Vectors::setValidity2

## From Biostrings ---------------------------------------
#' @importFrom Biostrings DNAStringSet
#' @export
Biostrings::DNAStringSet

#' @importFrom Biostrings AAStringSet
#' @export
Biostrings::AAStringSet

#' @importFrom Biostrings readDNAMultipleAlignment
#' @export
Biostrings::readDNAMultipleAlignment

#' @importFrom Biostrings translate
#' @export
Biostrings::translate

#' @importFrom Biostrings GENETIC_CODE_TABLE
#' @export
Biostrings::GENETIC_CODE_TABLE

#' @importFrom Biostrings getGeneticCode
#' @export
Biostrings::getGeneticCode


## From BiocGenerics ---------------------------------------
#' @importFrom BiocGenerics width
#' @export
BiocGenerics::width

#' @importFrom BiocGenerics start
#' @export
BiocGenerics::start

#' @importFrom BiocGenerics start<-
#' @export
BiocGenerics::`start<-`

#' @importFrom BiocGenerics end
#' @export
#' @examples 
#' ## Load an Automorphism object and get some 'end' coordinates
#' data("autm", package = "GenomAutomorphism")
#' end(autm[20:50])
BiocGenerics::end

#' @importFrom BiocGenerics end<-
#' @export
BiocGenerics::`end<-`

#' @importFrom BiocGenerics strand
#' @export
BiocGenerics::strand

#' @importFrom BiocGenerics strand<-
#' @export
BiocGenerics::`strand<-`


## From GenomicRanges ---------------------------------------
#' @importFrom GenomicRanges GRangesList
#' @export
GenomicRanges::GRangesList

## From numbers ---------------------------------------

#' @importFrom numbers modq
#' @export
numbers::modq

#' @importFrom numbers modlin
#' @export
numbers::modlin

