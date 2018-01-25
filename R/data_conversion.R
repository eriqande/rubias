

#' Convert Two-Column Genetic Data to Long Format
#'
#' Takes a data frame consisting of metadata followed by paired columns of genetic data,
#' with each column in the pair representing a gene copy at a locus.
#' Returns a list of two data frames: one with genetic data condensed into one column,
#' and the other with two-column structure intact, but with cleaned allele names.
#'
#' @keywords internal
#'
#' @param D A data frame containing two-column genetic data, optionally preceded by metadata.
#' The header of the first genetic data column in each pair lists the locus name, the second is ignored.
#' \strong{Locus names must not have spaces in them!}
#' @param gen_start_col The index (number) of the column in which genetic data starts.
#' Columns must be only genetic data after genetic data starts.
#'
#' @return \code{tcf2long} returns a list of two data frames: in the first, "long", the rightmost
#' column is the genetic data. Two new columns, "locus" and "gene copy", duplicate the original
#' column name provided in the first of each pair, and designate copies "a" and "b", respectively.
#' Metadata is duplicated as necessary for each locus. The second, "clean_short", replicates the
#' input dataset, but with column names replaced by "(locus name) a" and "(locus name) b" in each pair.
#' In other words the locus name has an "a" or a "b" added to it \emph{after a space}.
#'
#'
#' @examples
#' ## Convert the alewife dataset for further processing
#' ale_long <- tcf2long(alewife,17)
#' @export
tcf2long <- function(D, gen_start_col) {
  n = ncol(D)
  loc1 = seq(gen_start_col, n, by = 2)
  loc2 = seq(gen_start_col + 1, n, by = 2)
  dnam <- names(D)

  # checking requirements for gather: even number of loci, no duplicated or empty names, no spaces in names
  if( (n - gen_start_col) %% 2 == 0) stop("Odd number of locus columns. Check gen_start_col")
  nt <- table(dnam[loc1])
  nt <- nt[nt > 1]
  if(length(nt)>0) stop("Duplicated locus names: ", paste(names(nt), collapse = ", "))
  if(any(dnam=="") || any(is.na(dnam))) stop("Unnamed locus in genetic data")
  ns <- stringr::str_detect(dnam[loc1], " " )
  if(length(nt)>0) stop("Duplicated locus names: ", paste(names(nt), collapse = ", "))
  if( sum(ns) > 0) stop("Locus names cannot contain spaces: ",paste("\"", dnam[loc1][ns], "\"", sep = "", collapse = ", "))

  # here is a different way to do it that should be a little bit faster.
  num_loci <- length(loc1)
  num_short <- nrow(D)
  long <- tibble::tibble(
    sample_type = rep(D$sample_type, times = num_loci * 2),
    repunit = rep(D$repunit, times = num_loci * 2),
    collection = rep(D$collection, times = num_loci * 2),
    indiv = rep(D$indiv, times = num_loci * 2),
    locus = rep(dnam[loc1], each = num_short * 2),
    gene_copy = rep(c(rep("a", num_short), rep("b", num_short)), times = num_loci),
    allele = as.vector(as.matrix(D[, -(1:(gen_start_col - 1))]))
  )


  # renaming genetic column data
  names(D)[loc2] <- paste(names(D)[loc1], "b", sep = " ")
  names(D)[loc1] <- paste(names(D)[loc1], "a", sep = " ")

# long <- tidyr::gather_(data = D, key_col = "Locus", value_col = "allele", gather_cols = names(D)[gen_start_col:n]) %>%
#  tidyr::separate(col = Locus, into = c("locus","gene_copy"), sep = " ")



  list(long = long, clean_short = D)
}

#' Tabulate occurrences of all observed alleles in reference genetic data
#'
#' Takes the first output of \code{tcf2long}, along with two columns named "collection" and "sample_type",
#' and returns a data frame of allele counts for each locus within each reference population.
#' Alleles to be counted are identified from both reference and mixture populations.
#'
#' The "collection" column should be a key assigning samples to the desired groups,
#' e.g. collection site, run time, year.
#' The "sample_type" column must contain either "reference" or "mixture" for each sample.
#'
#' @keywords internal
#'
#' @param D A data frame containing, at minimum, a column of sample group identifiers named
#' "collection", a column designating each row as "reference" or "mixture", named "sample_type",
#' and (from tcf2long output) locus, gene copy, and observed alleles. If higher-level reporting
#' unit counts are desired, must have a column of reporting unit identifiers named "repunit"
#' @param pop_level a character vector expressing the population level for which allele counts
#' should be tabulated. Set to "collection" for collection/underlying sample group (default),
#' or "repunit" for reporting unit/overlying sample groups
#'
#' @return \code{reference_allele_counts} returns a long-format dataframe, with count data for
#' each collection, locus, and allele. Counts are only drawn from "reference" samples; alleles
#' unique to the "mixture" samples will still appear in the list, but will have 0s for all groups.
#'
#' @examples
#' ## count alleles in alewife reference populations
#' example(tcf2long)  # gets variable ale_long
#' ale_rac <- reference_allele_counts(ale_long$long)
#'
#'@export
reference_allele_counts <- function(D, pop_level = "collection") {
  # gathering the names of the gene loci and alleles across
  # all reference and mixture collections
  if(pop_level == "collection") {
    if( any(names(D) == "collection") == F) stop("no 'collection' column given")
  } else if(pop_level == "repunit") {
    if( any(names(D) == "repunit") == F) stop("no 'repunit' column given")
  } else stop("invalid selection for pop_level")
  if( any(names(D) == "sample_type") == F) stop("no 'sample_type' column given")
  if( any(is.na(D$collection)) ) stop("collection column may not contain NAs")
  if( any(is.na(D$sample_type)) ) stop("Sample_type column may not contain NAs")
  if(length(setdiff( unique(D$sample_type), c("reference","mixture"))) != 0) {
  warning("sample_type should not hold values other than 'reference' and 'mixture';
          other values ignored in allele counting")
  }

  D$collection <- factor(D$collection, levels = unique(D$collection))
  D$repunit <- factor(D$repunit, levels = unique(D$repunit))
  # generate data frame of alleles for each locus
  loc_alle_names <- D %>%
    dplyr::group_by(locus, allele) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(allele)) %>%
    dplyr::select(-n)

  if(pop_level == "collection") {
    # gathering the observed counts of each allele in each reference population
    ref_col_alle <- D %>%
      dplyr::filter(sample_type == "reference") %>%
      dplyr::group_by(collection, locus, allele) %>%
      dplyr::tally() %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(allele)) %>%
      dplyr::rename(counts = n)

    # creating combined table of all alleles & loci in each population,
    # whether in reference samples or not
    lefties <- lapply(unique(ref_col_alle$collection), function(x) cbind(collection = x, loc_alle_names, stringsAsFactors = FALSE)) %>%
      dplyr::bind_rows()

    # counts of each allele pair
    pop_tab <- suppressMessages(dplyr::left_join(lefties, ref_col_alle)) %>%
      dplyr::select(collection, dplyr::everything()) %>%
      dplyr::arrange(collection, locus, allele)

    pop_tab$counts[is.na(pop_tab$counts)] <- 0
    pop_tab
  } else {
    # gathering the observed counts of each allele in each reference population
    ref_col_alle <- D %>%
      dplyr::filter(sample_type == "reference") %>%
      dplyr::group_by(repunit, locus, allele) %>%
      dplyr::tally() %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(allele)) %>%
      dplyr::rename(counts = n)

    # creating combined table of all alleles & loci in each population,
    # whether in reference samples or not
    lefties <- lapply(unique(ref_col_alle$repunit), function(x) cbind(repunit = x, loc_alle_names, stringsAsFactors = FALSE)) %>%
      dplyr::bind_rows()

    # counts of each allele paired
    pop_tab <- suppressMessages(dplyr::left_join(lefties, ref_col_alle)) %>%
      dplyr::select(repunit, dplyr::everything()) %>%
      dplyr::arrange(repunit, locus, allele)

    pop_tab$counts[is.na(pop_tab$counts)] <- 0
    pop_tab
  }


}


#' Convert data frame of allele frequencies to nested lists
#'
#' List-izes the output of \code{reference_allele_counts} into a usable format for \code{allelic_list}
#'
#' @keywords internal
#'
#' @param D the long-format dataframe of counts by collection, locus, and allele,
#' output by \code{reference_allele_counts}, to be made into a nested list
#'
#' @return \code{a_freq_list} returns a list named by loci, each element of which is a matrix
#' containing that locus's allele count data. Rows in the matrix mark alleles, and columns collections
#'
#' @examples
#'
#'  # Generate a list of individual genotypes by allele from
#'  # the alewife data's reference allele count tables
#'  example(reference_allele_counts)
#'  ale_ac <- a_freq_list(ale_rac)
#' @export
a_freq_list <- function(D, pop_level = "collection") {
  if(pop_level == "collection"){
    tmp <- D %>%
      dplyr::arrange(locus,collection,allele)

    split(tmp, f = tmp$locus) %>%         # creates the list
      lapply(function(x) {         # names the rows by allele, and each column by collection
        nalle <- length(unique(x$allele))
        ret <- matrix(x$counts, nrow = nalle)
        dimnames(ret) <- list(allele = unique(x$allele), pop = unique(x$collection))
        ret
      })
  } else {
    tmp <- D %>%
      dplyr::arrange(locus,repunit,allele)

    split(tmp, f = tmp$locus) %>%         # creates the list
      lapply(function(x) {         # names the rows by allele, and each column by collection
        nalle <- length(unique(x$allele))
        ret <- matrix(x$counts, nrow = nalle)
        dimnames(ret) <- list(allele = unique(x$allele), pop = unique(x$repunit))
        ret
      })
  }

}


#' Create genotype lists for each locus
#'
#' Uses the allele counts from \code{a_freq_list} and the cleaned short-format output of
#' \code{tcf2long} to generate a nested list of individual genotypes for each locus
#'
#' @keywords internal
#'
#' @param cs a clean short genetic data matrix; the second element of the
#' output from \code{tcf2long}. Must have a column of individual identifiers, named "indiv"
#' @param ac allele counts from a_freq_list
#' @param samp_type choose which sample types of individuals to include in output:
#' "mixture", "both", or "reference"
#'
#' @return \code{allelic_list} returns a two-component nested list, with data stored as character
#' names of alleles ($chr) or as integer indices for the alleles ($int). Both forms contain lists
#' representing to loci, with two component vectors corresponding to gene copies a and b.
#'
#' @examples
#' example(a_freq_list)
#' ale_cs <- tcf2long(alewife, 17)$clean_short
#' # Get the vectors of gene copies a and b for all loci in integer index form
#' ale_alle_list <- allelic_list(ale_cs, ale_ac)$int
#'
#' @export
allelic_list <- function(cs, ac, samp_type = "both") {

  #Choosing a sample type causes the data set to be filtered by this type during allele list retrieval
  if(samp_type != "both") {
    cs <- cs %>%
      dplyr::filter(sample_type == samp_type)
  }
  # first, pull the cleaned first alleles for each locus into a list of character vectors
  all_a <- lapply(paste(names(ac), "a"), function(x) {
    ret <- as.character(cs[[x]])
    names(ret) <- cs$indiv
    ret
  })
  names(all_a) <- names(ac)

  # repeat with second gene copy
  all_b <- lapply(paste(names(ac), "b"), function(x) {
    ret <- as.character(cs[[x]])
    names(ret) <- cs$indiv
    ret
  })
  names(all_b) <- names(ac)

  # combine into list of locus lists, each of 2 character vectors, a and b
  ret <- lapply(names(ac), function(n) list(a = all_a[[n]], b = all_b[[n]]))
  names(ret) <- names(ac)

  #convert to factor, then to integers, while preserving sample IDs as names
  iret <- lapply(names(ac), function(n) {
    fa <- ret[[n]]$a
    fb <- ret[[n]]$b
    levs <- rownames(ac[[n]])
    a <- as.integer(factor(fa, levels = levs))
    b <- as.integer(factor(fb, levels = levs))
    names(a) <- names(fa)
    names(b) <- names(fb)
    list(a = a, b = b)
  })
  names(iret) <- names(ac)

  list(chr = ret, int = iret)
}


#' Collect essential data values before mixture proportion estimation
#'
#' Takes all relevant information created in previous steps of data conversion pipeline,
#' and combines into a single list which serves as input for further calculations
#'
#' Genotypes represented in \code{I_list} are converted into a single long vector,
#' ordered by locus, individual, and gene copy, with \code{NA} values represented as 0s.
#' Similarly, \code{AC_list} is unlisted to \code{AC}, ordered by locus, collection,
#' and allele. \code{DP} is a list of Dirichlet priors for likelihood calculations, created
#' by adding the values calculated from \code{alle_freq_prior} to each allele
#' \code{sum_AC} and \code{sum_DP} are the summed allele values for each locus
#' of their parent vectors, ordered by locus and collection.
#'
#' @keywords internal
#'
#' @param AC_list a list of allele count matrices; output from \code{a_freq_list}
#' @param I_list a list of genotype vectors; output from \code{allelic_list}
#' @param PO a vector of collection (population of origin) indices
#' for every individual in the sample, in order identical to \code{I_list}
#' @param coll_N a vector of the total number of individuals in each collection,
#' in order of appearance in the dataset
#' @param RU_vec a vector of collection indices, sorted by reporting unit
#' @param RU_starts a vector of indices, designating the first collection for each
#' reporting unit in RU_vec
#' @param alle_freq_prior a one-element named list specifying the prior to be used when
#' generating Dirichlet parameters for genotype likelihood calculations. The name of the
#' list item determines the type of prior used, with options \code{"const"}, \code{"scaled_const"},
#' and \code{"empirical"}. If \code{"const"}, the listed number will be taken as a constant
#' added to the count for each allele, locus, and collection.
#' If \code{"scaled_const"}, the listed number will be divided by the number of alleles at a locus,
#' then added to the allele counts. If \code{"empirical"}, the listed number will be multiplied
#' by the relative frequency of each allele across all populations, then added to the allele counts.
#'
#' @return \code{list_diploid_params} returns a list of the information necessary
#' for the calculation of genotype likelihoods in MCMC:
#'
#' \code{L}, \code{N}, and \code{C} represent the number of loci, individual genotypes,
#' and collections, respectively. \code{A} is a vector of the number of alleles at each
#' locus, and \code{CA} is the cumulative sum of \code{A}. \code{coll}, \code{coll_N},
#' \code{RU_vec}, and \code{RU_starts} are copied directly from input.
#'
#' \code{I}, \code{AC}, \code{sum_AC}, \code{DP}, and \code{sum_DP} are vectorized
#' versions of data previously represented as lists and matrices; indexing macros
#' use \code{L}, \code{N}, \code{C}, \code{A}, and \code{CA} to access these vectors
#' in later Rcpp-based calculations.
#'
#' @examples
#' \dontrun{
#' example(allelic_list)
#' PO <- as.integer(factor(ale_long$clean_short$collection))
#' coll_N <- as.vector(table(PO))
#'
#' Colls_by_RU <- dplyr::count(ale_long$clean_short, repunit, collection) %>%
#'    dplyr::select(-n)
#'  PC <- rep(0, length(unique((Colls_by_RU$repunit))))
#'  for(i in 1:nrow(Colls_by_RU)) {
#'    PC[Colls_by_RU$repunit[i]] <- PC[Colls_by_RU$repunit[i]] + 1
#'  }
#' RU_starts <- c(0, cumsum(PC))
#' RU_vec <- as.integer(Colls_by_RU$collection)
#' param_list <- list_diploid_params(ale_ac, ale_alle_list, PO, coll_N, RU_vec, RU_starts)
#' }
#'
#' @export
list_diploid_params <- function(AC_list, I_list, PO, coll_N, RU_vec, RU_starts,
                                alle_freq_prior = list("const_scaled" = 1)) {

  #check for a valid option in computing allele frequency priors
  if (!(names(alle_freq_prior)[1] %in% c("const_scaled", "const", "empirical"))) stop("Choice ", names(alle_freq_prior)[1], " unknown for allele frequency prior.")

  # get number of alleles at each locus
  as <- lapply(AC_list, nrow)

  loc_cycle <- names(AC_list)
  if(names(alle_freq_prior[1]) == "const_scaled") {
    DP_list <- lapply(loc_cycle, function(x) AC_list[[x]] + alle_freq_prior[[1]]/as[[x]])
  }
  if(names(alle_freq_prior[1]) == "const") {
    DP_list <- lapply(loc_cycle, function(x) AC_list[[x]] + alle_freq_prior[[1]])
  }
  if(names(alle_freq_prior[1]) == "empirical") {
    DP_list <- lapply(AC_list, function(x){
      rel_freq <- rowSums(x)/sum(x)
      x + (rel_freq * alle_freq_prior[[1]])
    })
  }

  list(L = length(AC_list),
       N = length(I_list[[1]]$a),
       C = ncol(AC_list[[1]]),
       A = unname(sapply(AC_list, nrow)),
       CA = unname(c(0, cumsum(sapply(AC_list, nrow)))),
       coll = PO,
       coll_N = coll_N,
       RU_vec = RU_vec,
       RU_starts = RU_starts,
       I = unname(unlist(lapply(I_list, function(l) {
         x <- t(cbind(l$a, l$b))
         x[is.na(x)] <- 0
         x
       }))),
       AC = unname(unlist(AC_list)),
       sum_AC = lapply(AC_list, colSums) %>%
         unlist(),
       DP = unlist(DP_list),
       sum_DP = lapply(DP_list,colSums) %>%
         unlist())
}

#' Generate MCMC parameter list from two-column genetic data & print summary
#'
#' This function is a wrapper for all steps to create the parameter list necessary
#' for genotype log-likelihood calculation from the starting two-column genetic data
#'
#' In order for all steps in conversion to be carried out successfully, the dataset
#' must have "repunit", "collection", "indiv", and "sample_type" columns preceding
#' two-column genetic data. If \code{summ == TRUE}, the function prints summary statistics
#' describing the structure of the dataset, as well as the presence of missing data,
#' enabling verification of proper data conversion.
#'
#' @param D A data frame containing two-column genetic data, preceded by metadata.
#' The header of the first genetic data column in each pair lists the locus name,
#' the second is ignored. \strong{Locus names must not have spaces in them!}
#' Required metadata includes a column of unique individual identifiers named "indiv",
#' a column named "collection" designating the sample groups, a column "repunit"
#' designating the reporting unit of origin of each fish, and a "sample_type" column
#' denoting each individual as a "reference" or "mixture" sample. \emph{No NAs should be
#' present in metadata}
#' @param gen_start_col The index (number) of the column in which genetic data starts.
#' Columns must be only genetic data after genetic data starts.
#' @param samp_type the sample groups to be include in the individual genotype list,
#' whose likelihoods will be used in MCMC. Options "reference", "mixture", and "both"
#' @param alle_freq_prior a one-element named list specifying the prior to be used when
#' generating Dirichlet parameters for genotype likelihood calculations. Valid methods
#' include \code{"const"}, \code{"scaled_const"}, and \code{"empirical"}. See
#' \code{?list_diploid_params} for method details.
#' @param summ logical indicating whether summary descriptions of the formatted data be provided
#'
#' @return \code{tcf2param_list} returns the output of \code{list_diploid_params},
#' after the original dataset is converted to a usable format and all relevant values
#' are extracted. See \code{?list_diploid_params} for details
#'
#' @keywords internal
#' @examples
#' ale_par_list <- tcf2param_list(alewife, 17)
#'
#' @export
tcf2param_list <- function(D, gen_start_col, samp_type = "both",
                           alle_freq_prior = list("const_scaled" = 1), summ = T){

  # coerce collection and repunit to be factors.  This is important since things get turned
  # into integers
  D$collection <- factor(D$collection, levels = unique(D$collection))
  D$repunit <- factor(D$repunit, levels = unique(D$repunit))

  cleaned <- tcf2long(D, gen_start_col)
  AC_list <- reference_allele_counts(cleaned$long) %>%
    a_freq_list()
  I_list <- allelic_list(cleaned$clean_short, AC_list, samp_type = samp_type)$int
  #PO <- as.integer(factor(cleaned$clean_short$collection))
  if(samp_type == "both") {
    PO <- as.integer(cleaned$clean_short$collection)
#    PO <- dplyr::select(cleaned$clean_short, collection) %>%
#      simplify2array() %>%
#      factor(., levels = unique(D$collection)) %>%
#      as.integer()

    coll_N <- unname(table(cleaned$clean_short$collection))
#   coll_N <- dplyr::count(cleaned$clean_short, collection) %>%
#      dplyr::select(n)
#    coll_N <- unname(unlist(coll_N))
  } else {
    PO <- dplyr::filter(cleaned$clean_short, sample_type == samp_type) %>%
      dplyr::select(collection) %>%
      simplify2array() %>%
      factor(., levels = unique(D$collection)) %>%
      as.integer()
    coll_N <- dplyr::filter(cleaned$clean_short, sample_type == samp_type) %>%
      dplyr::count(collection) %>%
      dplyr::select(n)
    coll_N <- unname(unlist(coll_N))
  }

  ### ONLY designed for factorized repunit and collection
  Colls_by_RU <- dplyr::count(cleaned$clean_short[, c("repunit", "collection")], repunit, collection) %>%
    dplyr::select(-n)
#  Colls_by_RU <- dplyr::count(cleaned$clean_short, repunit, collection) %>%
#    dplyr::select(-n)
  PC <- rep(0, length(unique((Colls_by_RU$repunit))))
  for(i in 1:nrow(Colls_by_RU)) {
    PC[Colls_by_RU$repunit[i]] <- PC[Colls_by_RU$repunit[i]] + 1
  }
  RU_starts <- c(0, cumsum(PC))
  names(RU_starts) <- as.character(unique(Colls_by_RU$repunit))
  RU_vec <- as.integer(Colls_by_RU$collection)
  names(RU_vec) <- as.character(Colls_by_RU$collection)
  params <- list_diploid_params(AC_list, I_list, PO, coll_N, RU_vec, RU_starts, alle_freq_prior)
  percent.missing <- sum(params$I == 0)/length(params$I) * 100
  RU_list <- unique(cleaned$clean_short$repunit)

  # here, the names of the indivs, collections, and repunits in the order in which they
  # appear as integers in this data structure.  We should have done this on day one, but
  # we didn't which is why there are so many mysterious ways that have been used to get the
  # names of the collections or repunits.  At some point, we will have to go through
  # and standardize all those different name accesses. But not now...
  params$indiv_names <- D$indiv
  params$collection_names <- levels(D$collection)
  params$repunit_names <- levels(D$repunit)
  params$locus_names <- names(AC_list)

  if(summ == T){
    cat(paste('Summary Statistics:',
              paste(params$N, 'Individuals in Sample'),
              paste(paste(params$L, 'Loci:'), paste(names(AC_list), collapse = ", ")),
              paste(paste(paste(length(RU_list), 'Reporting Units:'),
                          paste(levels(RU_list), collapse = ", "))),
              paste(paste(paste(params$C, 'Collections:'),
                          paste(colnames(AC_list[[1]]), collapse = ", "))),
              paste(percent.missing, '% of allelic data identified as missing', sep = ""),
              sep = '\n\n'))
  }

  params
}


#' Get the average within-RU assignment rate for each collection
#'
#' This function takes a matrix of scaled genotype likelihoods for a group of
#' individuals of known origin, and calculates the average rate at which individuals
#' in a particular collection are assigned to the correct reporting unit.
#'
#' The average rate of correct within-reporting unit assignment is proportional to
#' reporting-unit-level bias in the posterior probability for this collection;
#' if the correct assignment rate is high relative to other collections,
#' it will be upwardly biased, and vice versa. The inverse of this vector is used
#' to scale Dirichlet draws of \code{omega} during misassignment-scaled MCMC.
#' @param SL a scaled likelihood matrix; each column should sum to one, and represent
#' the probability of assignments to each collection (row) for a particular individual
#' @param coll a vector of the collection of origin indices of the individuals (length = \code{ncol(SL)})
#' @param RU_starts a vector delineating starting indices of different reporting units in RU_vec
#' @param RU_vec a vector of collection indices, organized by reporting unit
#'
#' @return \code{avg_coll2correctRU} returns a vector of length \code{nrow(SL)}, where
#' each element represents the average proportion of fish from the corresponding collection
#' which are correctly assigned to the proper collection, or misassigned to another
#' collection within the same reporting unit. This is distinct from the rate of
#' correct assignment at the collection level, which is too low and variable to serve
#'  as a stable metric for \code{omega} scaling.
#' @examples
#' params <- tcf2param_list(alewife, 17)
#' SL <- geno_logL(params) %>% exp() %>% apply(2, function(x) x/sum(x))
#' correct <- avg_coll2correctRU(SL, params$coll, params$RU_starts, params$RU_vec)
#' @export
#' @keywords internal
avg_coll2correctRU <- function(SL, coll, RU_starts, RU_vec) {
  ru_SL <- lapply(1:(length(RU_starts)-1), function(ru){
    apply(rbind(SL[RU_vec[(RU_starts[ru]+1):RU_starts[ru + 1]], ]), 2, sum)
  }) %>%
    do.call("rbind", .)
  avg_SL <- lapply(1:nrow(SL), function(x) as.array(ru_SL[, (coll == x)])) %>%
    lapply(function(x) apply(x,MARGIN = 1, FUN = mean)) %>%
    simplify2array()
  correct <- lapply(1:(length(RU_starts)-1), function(ru) {
    ru_colls <- RU_vec[(RU_starts[ru]+1):RU_starts[ru+1]]
    out <- avg_SL[ru, ru_colls]
  }) %>%
    unlist()
  correct
}
