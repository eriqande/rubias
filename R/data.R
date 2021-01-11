
#' Microsat data from alewife herring reference populations
#'
#' Standard two-column genetic data with lots of other columns preceding it.
#' Can be fed directly into rubias because it has at least the columns sample_type,
#' collection, repunit and indiv.
#' @format A tibble.
#' @source \url{https://datadryad.org/stash/dataset/doi:10.5061/dryad.80f4f}
#' @docType data
#' @name alewife
NULL


#' Microsat data from blueback herring reference populations
#'
#' Standard two-column genetic data with lots of other columns preceding it.
#' Can be fed directly into rubias because it has at least the columns sample_type,
#' collection, repunit and indiv.
#' @format A tibble.
#' @source \url{https://datadryad.org/stash/dataset/doi:10.5061/dryad.80f4f}
#' @docType data
#' @name blueback
NULL




#' SNP data from chinook reference populations
#'
#' Chinook salmon baseline data similar to that which can be
#' downloaded from \url{https://datadryad.org/stash/dataset/doi:10.5061/dryad.574sv}.
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
#' @source \url{https://datadryad.org/stash/dataset/doi:10.5061/dryad.574sv}
#' @docType data
#' @name chinook
NULL


#' SNP data from selected chinook reference populations
#'
#' A small number of poulations from the
#' Chinook salmon baseline data similar to that which can be
#' downloaded from \url{https://datadryad.org/stash/dataset/doi:10.5061/dryad.574sv}.
#' This data set includes 91 SNPs and 909 fish.
#' @format A tbl_df-ed (from dplyr) data frame with 909 rows and 185 variables. The first
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
#' @source \url{https://datadryad.org/stash/dataset/doi:10.5061/dryad.574sv}
#' @docType data
#' @name small_chinook_ref
NULL


#' SNP data from Chinook salmon taken in May/August 2015 from California fisheries
#'
#' This has data from 91 SNP markers (a subset of the 95 markers in the \code{\link{chinook}} baseline
#' data set).
#' @format A tbl_df-ed (from dplyr) data frame with 2256 rows and 193 variables. The first
#' four columns are meta data.
#' The remaining columns are two columns for each locus.  These columns are named like,
#' "Locus.1" and "Locus.2" for the first and second gene copies at that locus.  For example,
#' "Ots_104569-86.1" and  "Ots_104569-86.2".  The locus columns are ints and missing data
#' is denoted by NA.
#' @source Southwest Fisheries Science Center, Santa Cruz, CA
#' @docType data
#' @name chinook_mix
NULL


#' Small sample of SNP data from Chinook salmon taken in May/August 2015 from California fisheries
#'
#' This is simply a sample of 100 fish from \code{\link{chinook}}.
#'
#' @source Southwest Fisheries Science Center, Santa Cruz, CA
#' @docType data
#' @name small_chinook_mix
NULL



#' List of example ways of specifying repunit and collection quantities in simulations
#'
#' This is just a list of tibbles that can be passed to the alpha_repunit or the
#' alpha_collection parameters in, for example, \code{\link{assess_reference_loo}}.
#' @source  Made it up!
#' @docType data
#' @name sim_spec_examples
NULL


#' perfect-assignment genetic data for chinook.
#'
#' This is just like the \code{\link{chinook}} data, but only has 7 loci and all loci are
#' fixed in fortuitous patterns so that every single collection is easily resolved.  This is
#' primarily useful for testing purposes.
#' @source Made it up!
#' @docType data
#' @name perfect_chinook
NULL


#' perfect-assignment mixture genetic data for chinook.
#'
#' This is similar to the \code{\link{chinook_mix}} data, but only has 7 loci and all loci are
#' fixed in fortuitous patterns so that every single collection is easily resolved.  This is
#' primarily useful for testing purposes.  The name of the individual has its collection
#' inside the colons.
#' @source Made it up!
#' @docType data
#' @name perfect_chinook_mix
NULL


#' a vector that gives a desired sort order of the chinook collections
#'
#' This is just an example of what one would use as levels in order to
#' get the \code{\link{chinook}} collections in a desired sort order after
#' analysis.  The issue here is collection in the input data frame to most
#' functions must be a character vector, not a factor.  But, after analysis
#' you can always make them a factor again and use a vector like this
#' one to specify the levels.
#' @source Made it up!
#' @docType data
#' @name chinook_collection_levels
NULL



#' a vector that gives a desired sort order of the chinook repunits
#'
#' This is just an example of what one would use as levels in order to
#' get the \code{\link{chinook}} repunits in a desired sort order after
#' analysis.  The issue here is that repunit in the input data frame to most
#' functions must be a character vector, not a factor.  But, after analysis
#' you can always make them a factor again and use a vector like this
#' one to specify the levels.
#' @source Made it up!
#' @docType data
#' @name chinook_repunit_levels
NULL

