
#' Microsat data from alewife herring reference populations
#'
#' gotta say more about it.  And hey! We need to put the reporting unit designations on it.
#' @format A tbl_df-ed (from dplyr) data frame with XXX rows and YY variables. The relevant variables for
#' analyses here are:
#' \describe{
#'   \item{Field_Num}{unique identifier for each bird}
#'   \item{Isotope.Value}{hydrogen isotope ratios measured in the bird's feather}
#'   \item{lat}{latitude of the bird's breeding/sampling location}
#'   \item{long}{latitude of the bird's breeding/sampling location}
#' }
#' @source Got it off GitHub, but must be at Dryad now too...
#' @docType data
#' @name alewife
NULL


#' Microsat data from blueback herring reference populations.
#'
#' gotta say more about it
#' @format A tbl_df-ed (from dplyr) data frame with XXX rows and YY variables. The relevant variables for
#' analyses here are:
#' \describe{
#'   \item{Field_Num}{unique identifier for each bird}
#'   \item{Isotope.Value}{hydrogen isotope ratios measured in the bird's feather}
#'   \item{lat}{latitude of the bird's breeding/sampling location}
#'   \item{long}{latitude of the bird's breeding/sampling location}
#' }
#' @source Got it off GitHub, but must be at Dryad now too...
#' @docType data
#' @name blueback
NULL




#' SNP data from chinook reference populations
#'
#' Chinook salmon baseline data similar to that which can be
#' downloaded from \link{http://datadryad.org/resource/doi:10.5061/dryad.574sv/1}.
#' This data set includes 91 SNPs and 7301 fish and is what the Dryad data became
#' after we converted from TaqMan to SNPtype assays (being forced to toss some loci)
#' and tossed out a bunch of lousy historical
#' samples from Trinity River.
#' @format A tbl_df-ed (from dplyr) data frame with 7,301 rows and 185 variables. The first
#' three columns are
#' \describe{
#'   \item{repunit (chr)}{the reporting unit that the individual is in}
#'   \item{pop (chr)}{the population from which the individual was sampled}
#'   \item{ID (chr)}{Unique identifier of the individual fish}
#' }
#' The remaining columns are two columns for each locus.  These columns are named like,
#' "Locus.1" and "Locus.2" for the first and second gene copies at that locus.  For example,
#' "Ots_104569-86.1" and  "Ots_104569-86.2".  The locus columns are ints and missing data
#' is denoted by NA.
#' @source \link{http://datadryad.org/resource/doi:10.5061/dryad.574sv/1}
#' @docType data
#' @name chinook
NULL


#' SNP data from Chinook salmon taken in May/August 2015 from California fisheries
#'
#' This has data from 91 SNP markers (a subset of the 95 markers in the \code{\link{chinook}} baseline
#' data set).
#' @format A tbl_df-ed (from dplyr) data frame with 8,031 rows and 193 variables. The first
#' three columns are
#' \describe{
#'   \item{ID (chr)}{Unique identifier of the individual fish}
#' }
#' The remaining columns are two columns for each locus.  These columns are named like,
#' "Locus.1" and "Locus.2" for the first and second gene copies at that locus.  For example,
#' "Ots_104569-86.1" and  "Ots_104569-86.2".  The locus columns are ints and missing data
#' is denoted by NA.
#' @source Southwest Fisheries Science Center, Santa Cruz, CA
#' @docType data
#' @name chinook_mix
NULL


#' SNP data from Chinook reference and mixtures mixed together
#'
#' gotta say more
#' @source Southwest Fisheries Science Center, Santa Cruz, CA
#' @docType data
#' @name chinook_all
#' NULL

