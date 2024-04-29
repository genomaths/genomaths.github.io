#' Reexport useful functions to be available to users
## From S4Vectors ---------------------------------------

#' @importFrom S4Vectors mcols
#' @returns The same as in \code{\link[S4Vectors]{mcols}}.
#' @export
#' @examples 
#' ## Load an Automorphism object and take its metacolumns
#' data("autm", package = "GenomAutomorphism")
#' mcols(autm)
S4Vectors::mcols

#' @importFrom S4Vectors mcols<-
#' @returns The same as in \code{\link[S4Vectors]{mcols}}. 
#' @export
S4Vectors::`mcols<-`

#' @importFrom S4Vectors setValidity2
#' @returns The same as in \code{\link[S4Vectors]{setValidity2}}. 
#' @export
S4Vectors::setValidity2

## From Biostrings ---------------------------------------
#' @importFrom Biostrings DNAStringSet
#' @returns The same as in \code{\link[Biostrings]{DNAStringSet}}. 
#' @export
Biostrings::DNAStringSet

#' @importFrom Biostrings AAStringSet
#' @returns The same as in \code{\link[Biostrings]{AAStringSet}}. 
#' @export
Biostrings::AAStringSet

#' @importFrom Biostrings readDNAMultipleAlignment
#' @returns The same as in \code{\link[Biostrings]{readDNAMultipleAlignment}}. 
#' @export
Biostrings::readDNAMultipleAlignment

#' @importFrom Biostrings DNAMultipleAlignment
#' @returns The same as in \code{\link[Biostrings]{DNAMultipleAlignment}}. 
#' @export
Biostrings::DNAMultipleAlignment

#' @importFrom Biostrings AAMultipleAlignment
#' @returns The same as in \code{\link[Biostrings]{AAMultipleAlignment}}. 
#' @export
Biostrings::AAMultipleAlignment

#' @importFrom XVector subseq
#' @returns The same as in \code{\link[XVector]{subseq}}. 
#' @export
XVector::subseq

#' @importFrom Biostrings translate
#' @returns The same as in \code{\link[Biostrings]{translate}}. 
#' @export
Biostrings::translate

#' @importFrom Biostrings GENETIC_CODE_TABLE
#' @returns The same as in \code{\link[Biostrings]{GENETIC_CODE_TABLE}}. 
#' @export
Biostrings::GENETIC_CODE_TABLE

#' @importFrom Biostrings getGeneticCode
#' @returns The same as in \code{\link[Biostrings]{getGeneticCode}}. 
#' @export
Biostrings::getGeneticCode


#' @importFrom Biostrings unmasked
#' @returns The same as in \code{\link[Biostrings]{unmasked}}. 
#' @export
Biostrings::unmasked

## From BiocGenerics ---------------------------------------
#' @importFrom BiocGenerics width
#' @returns The same as in \code{\link[BiocGenerics]{width}}. 
#' @export
BiocGenerics::width

#' @importFrom BiocGenerics start
#' @returns The same as in \code{\link[BiocGenerics]{start}}. 
#' @export
BiocGenerics::start

#' @importFrom BiocGenerics start<-
#' @returns The same as in \code{\link[BiocGenerics]{start}}. 
#' @examples
#' ## See \code{\link[BiocGenerics]{start}}.
#' 
#' @export
BiocGenerics::`start<-`

#' @importFrom BiocGenerics end
#' @returns The same as in \code{\link[BiocGenerics]{end}}. 
#' @export
#' @examples 
#' ## Load an Automorphism object and get some 'end' coordinates
#' data("autm", package = "GenomAutomorphism")
#' end(autm[20:50])
BiocGenerics::end

#' @importFrom BiocGenerics end<-
#' @returns The same as in \code{\link[BiocGenerics]{end}}. 
#' @export
BiocGenerics::`end<-`

#' @importFrom BiocGenerics strand
#' @returns The same as in \code{\link[BiocGenerics]{strand}}. 
#' @export
BiocGenerics::strand

#' @importFrom BiocGenerics strand<-
#' @returns The same as in \code{\link[BiocGenerics]{strand}}. 
#' @export
BiocGenerics::`strand<-`


## From GenomicRanges ---------------------------------------
#' @importFrom GenomicRanges GRangesList
#' @returns The same as in \code{\link[GenomicRanges]{GRangesList}}. 
#' @export
GenomicRanges::GRangesList

#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @returns The same as in 
#' \code{\link[GenomicRanges]{makeGRangesFromDataFrame}}. 
#' @export
GenomicRanges::makeGRangesFromDataFrame

## From GenomicRanges ---------------------------------------

#' @importFrom numbers modq
#' @returns The same as in \code{\link[numbers]{modq}}. 
#' @export
numbers::modq

#' @importFrom numbers modlin
#' @returns The same as in \code{\link[numbers]{modlin}}. 
#' @export
numbers::modlin

## From matrixStats ---------------------------------------

#' @importFrom matrixStats rowSums2
#' @returns The same as in \code{\link[matrixStats]{rowSums2}}. 
#' @export
matrixStats::rowSums2

#' @importFrom matrixStats colSums2
#' @returns The same as in \code{\link[matrixStats]{colSums2}}. 
#' @export
matrixStats::colSums2

#' @importFrom matrixStats colMeans2
#' @returns The same as in \code{\link[matrixStats]{colMeans2}}. 
#' @export
matrixStats::colMeans2

#' @importFrom matrixStats rowMeans2
#' @returns The same as in \code{\link[matrixStats]{rowMeans2}}. 
#' @export
matrixStats::rowMeans2

#' @importFrom matrixStats rowVars
#' @returns The same as in \code{\link[matrixStats]{rowVars}}. 
#' @export
matrixStats::rowVars

#' @importFrom matrixStats colVars
#' @returns The same as in \code{\link[matrixStats]{colVars}}. 
#' @export
matrixStats::colVars

#' @importFrom matrixStats colSds
#' @returns The same as in \code{\link[matrixStats]{colSds}}. 
#' @export
matrixStats::colSds

#' @importFrom matrixStats rowSds
#' @returns The same as in \code{\link[matrixStats]{rowSds}}. 
#' @export
matrixStats::rowSds













