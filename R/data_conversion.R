

#' Convert Two-Column Genetic Data to Long Format
#'
#' Takes a data frame consisting of metadata followed by paired columns of genetic data,
#' with each column in the pair representing a gene copy at a locus.
#' Returns a list of two data frames: one with genetic data condensed into one column,
#' and the other with two-column structure intact, but with cleaned allele names.
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
#' ale_long <- tcf2long(alewife,15)
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

  # renaming genetic column data
  names(D)[loc2] <- paste(names(D)[loc1], "b", sep = " ")
  names(D)[loc1] <- paste(names(D)[loc1], "a", sep = " ")

  long <- tidyr::gather_(data = D, key_col = "Locus", value_col = "allele", gather_cols = names(D)[gen_start_col:n]) %>%
    tidyr::separate(col = Locus, into = c("locus","gene_copy"), sep = " ")

  list(long = long, clean_short = D)
}

#' Tabulate occurences of all observed alleles in reference genetic data
#'
#' Takes the first output of \code{tcf2long}, along with two columns named "collection" and "sample_type",
#' and returns a data frame of allele counts for each locus within each reference population.
#' Alleles to be counted are identified from both reference and mixture populations.
#'
#' The "collection" column should be a key assigning samples to the desired groups,
#' e.g. collection site, run time, year.
#' The "sample_type" column must contain either "reference" or "mixture" for each sample.
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

    # counts of each allele paired
    pop_tab <- suppressMessages(dplyr::left_join(lefties, ref_col_alle)) %>%
      dplyr::select(collection, everything()) %>%
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
      dplyr::select(repunit, everything()) %>%
      dplyr::arrange(repunit, locus, allele)

    pop_tab$counts[is.na(pop_tab$counts)] <- 0
    pop_tab
  }


}


#' Convert data frame of allele frequencies to nested lists
#'
#' List-izes the output of \code{reference_allele_counts} into a usable format for allelic_list
#'
#'@param D the long-format dataframe of counts by collection, locus, and allele,
#'output by reference_allele_counts, to be made into a nested list
#'
#'@return \code{a_freq_list} returns a list named by loci, each element of which is a matrix
#'containing that locus's count data. Rows in the matrix mark alleles, and columns collections
#'
#'@examples#'
#' # Generate a list of individual genotypes by allele from
#' # the alewife data's reference allele count tables
#' example(reference_allele_counts)
#' ale_ac <- a_freq_list(ale_rac)
#'
#'@export
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
#' @param cs a clean short genetic data matrix, the second element of the
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
#' ale_cs <- tcf2long(alewife, 15)$clean_short
#' # Get the vectors of gene copies a and b for all loci in integer index form
#' ale_alle_list <- allelic_list(ale_cs, ale_ac)$int
#'
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

#' Get summed Dirichlet parameters for all alleles and populations in a numeric vector
#'
#' @param DP the list of dirichlet parameter matrices
#'
#' @examples
#' example(allelic_list)
#' as <- lapply(ale_ac, nrow)
#' loc_cycle <- names(ale_ac)
#' names(loc_cycle) <- loc_cycle
#' ale_DP_list <- lapply(loc_cycle, function(x) ale_ac[[x]] + 1/as[[x]])
#' sum_DP <- DP_sum(ale_DP_list)
#'
#' @export
DP_sum <- function(DP_list) {
  sum_DP <- lapply(DP_list, colSums)
  sum_DP <- unlist(sum_DP)
}


#' Collect essential data values before mixture proportion estimation
#'
#' @param AC_list output from reference_allele_counts
#' @param I_list output from allelic_list
#' @param PO vector of collection (population of origin) indices
#' for every individual in the sample
#' @param coll_N a vector of the total number of individuals in each collection
#' @param RU_vec a vector of collection indexes, sorted by reporting unit
#' @param RU_starts a vector of indexes, designating the first collection for each
#' reporting unit in RU_vec
#'
#' @return \code{list_diploid_params} returns a list of the information necessary
#' for the calculation of genotype likelihoods in MCMC
#'
#' @examples
#' example(allelic_list)
#' PO <- as.integer(factor(ale_long$clean_short$collection))
#' coll_N <- dplyr::count(ale_long$clean_short, collection) %>%
#'   dplyr::select(n)
#'   coll_N <- unname(unlist(coll_N))
#' Colls_by_RU <- dplyr::count(ale_long$clean_short, repunit, collection) %>%
#'    dplyr::select(-n)
#'  PC <- rep(0, length(unique((Colls_by_RU$repunit))))
#'  for(i in 1:nrow(Colls_by_RU)) {
#'    PC[Colls_by_RU$repunit[i]] <- PC[Colls_by_RU$repunit[i]] + 1
#'  }
#' RU_starts <- c(0, cumsum(PC))
#' RU_vec <- as.integer(Colls_by_RU$collection)
#' param_list <- list_diploid_params(ale_ac, ale_alle_list, PO, coll_N, RU_vec, RU_starts)
#'
#' @export
list_diploid_params <- function(AC_list, I_list, PO, coll_N, RU_vec, RU_starts) {
  as <- lapply(AC_list, nrow)
  loc_cycle <- names(AC_list)
  DP_list <- lapply(loc_cycle, function(x) AC_list[[x]] + 1/as[[x]])
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
       sum_AC = DP_sum(AC_list),
       DP = unlist(DP_list),
       sum_DP = DP_sum(DP_list))
}

#' Generate MCMC parameter list from two-column genetic data & Print Summary
#'
#' This function is a wrapper for all steps to create the parameter list necessary
#' for genotype log-likelihood calculation from the starting two-column genetic data
#'
#' @param D A data frame containing two-column genetic data, preceded by metadata.
#' The header of the first genetic data column in each pair lists the locus name,
#' the second is ignored. \strong{Locus names must not have spaces in them!}
#' Required metadata includes a column of unique individual identifiers named "indiv",
#' a column named "collection" designating the sample groups, and a "sample_type" column
#' denoting each individual as a "reference" or "mixture" sample
#' @param gen_start_col The index (number) of the column in which genetic data starts.
#' Columns must be only genetic data after genetic data starts.
#' @param samp_type the sample groups to be include in the individual genotype list,
#' whose likelihoods will be used in MCMC. Options "reference", "mixture", and "both"
#' @param summ (LOGICAL) should summary descriptions of the formatted data be provided?
#'
#' @examples
#' ale_par_list <- tcf2param_list(alewife, 15)
#'
#' @export
tcf2param_list <- function(D, gen_start_col, samp_type = "both", summ = T){

  #confirm that collection and repunit are factors; important for later processing
  D$collection <- factor(D$collection, levels = unique(D$collection))
  D$repunit <- factor(D$repunit, levels = unique(D$repunit))

  cleaned <- tcf2long(D, gen_start_col)
  AC_list <- reference_allele_counts(cleaned$long) %>%
    a_freq_list()
  I_list <- allelic_list(cleaned$clean_short, AC_list, samp_type = samp_type)$int
  #PO <- as.integer(factor(cleaned$clean_short$collection))
  if(samp_type == "both") {
    PO <- dplyr::select(cleaned$clean_short, collection) %>%
      simplify2array() %>%
      factor(., levels = unique(D$collection)) %>%
      as.integer()
    coll_N <- dplyr::count(cleaned$clean_short, collection) %>%
      dplyr::select(n)
    coll_N <- unname(unlist(coll_N))
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
  Colls_by_RU <- dplyr::count(cleaned$clean_short, repunit, collection) %>%
    dplyr::select(-n)
  PC <- rep(0, length(unique((Colls_by_RU$repunit))))
  for(i in 1:nrow(Colls_by_RU)) {
    PC[Colls_by_RU$repunit[i]] <- PC[Colls_by_RU$repunit[i]] + 1
  }
  RU_starts <- c(0, cumsum(PC))
  names(RU_starts) <- as.character(unique(Colls_by_RU$repunit))
  RU_vec <- as.integer(Colls_by_RU$collection)
  names(RU_vec) <- as.character(Colls_by_RU$collection)
  params <- list_diploid_params(AC_list, I_list, PO, coll_N, RU_vec, RU_starts)
  percent.missing <- sum(params$I == 0)/length(params$I) * 100
  RU_list <- unique(cleaned$clean_short$repunit)

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

#' generate a matrix of the average percentage of a collection assigned to each collection
#'
#' @param SL a scaled likelihood matrix
#' @param coll a vector of the collections of origin of the columns in SL
#'
#' @return a matrix in which each column represents the collection of origin, and each row
#' the average scaled likelihood of assignment to that collection
#'
#' @examples
#' params <- tcf2param_list(alewife, 15)
#' SL <- geno_logL(params) %>% exp() %>% apply(2, function(x) x/sum(x))
#' avg_SL <- avg_SL_coll(SL, params$coll)
#'
#' @export
avg_SL_coll <- function(SL, coll) {
  avg_SL <- lapply(1:nrow(SL), function(x) as.array(SL[, (coll == x)])) %>%
    lapply(function(x) apply(x,MARGIN = 1, FUN = mean)) %>%
    simplify2array()
  avg_SL
}

#' generate a matrix of the average percentage of a reporting unit misassigned
#'
#' @param SL a scaled likelihood matrix
#' @param params the MCMC parameter list output by list_diploid_params and tcf2param_list
#'
#' @examples
#' params <- tcf2param_list(alewife, 15)
#' SL <- geno_logL(params) %>% exp() %>% apply(2, function(x) x/sum(x))
#' miss <- avg_miss_ru(SL, params$coll, params$RU_starts, params$RU_vec)
#'
#' @export
avg_miss_ru <- function(SL, coll, RU_starts, RU_vec) {
  avg_SL <- lapply(1:nrow(SL), function(x) as.array(SL[, (coll == x)])) %>%
    lapply(function(x) apply(x,MARGIN = 1, FUN = mean)) %>%
    simplify2array()
  miss <- lapply(1:(length(RU_starts)-1), function(ru){
    n_coll <- RU_starts[ru+1] - RU_starts[ru]
    ru_colls <- RU_vec[(RU_starts[ru]+1):RU_starts[ru+1]]
    miss <- sum(avg_SL[-ru_colls, ru_colls])/n_coll
  }) %>%
    unlist
  names(miss) <- names(RU_starts[1:(length(RU_starts) - 1)])
  miss
}

#' generate a matrix of the average percentage of a collection misassigned to a different reporting unit
#'
#' @param SL a scaled likelihood matrix
#' @param coll a vector of the collection of origin indices of the fish (length == ncol(SL))
#' @param RU_starts a vector delineating starting indices of different reporting units in RU_vec
#' @param RU_vec a vector of collection indices, organized by reporting unit
#'
#' @examples
#' params <- tcf2param_list(alewife, 15)
#' SL <- geno_logL(params) %>% exp() %>% apply(2, function(x) x/sum(x))
#' miss <- avg_miss_ru(SL, params$coll, params$RU_starts, params$RU_vec)
#'
#' @export
avg_coll_miss_ru <- function(SL, coll, RU_starts, RU_vec) {
  ru_SL <- lapply(1:(length(RU_starts)-1), function(ru){
    apply(SL[RU_vec[(RU_starts[ru]+1):RU_starts[ru + 1]], ], 2, sum)
  }) %>%
    do.call("rbind", .)
  avg_SL <- lapply(1:nrow(SL), function(x) as.array(ru_SL[, (coll == x)])) %>%
    lapply(function(x) apply(x,MARGIN = 1, FUN = mean)) %>%
    simplify2array()
  miss <- lapply(1:(length(RU_starts)-1), function(ru) {
    ru_colls <- RU_vec[(RU_starts[ru]+1):RU_starts[ru+1]]
    miss <- colSums(avg_SL[-ru, ru_colls])
  }) %>%
    unlist()
  miss
}


#' generate a matrix of the average percentage of a collection assigned to the correct reporting unit
#'
#' @param SL a scaled likelihood matrix
#' @param coll a vector of the collection of origin indices of the fish (length == ncol(SL))
#' @param RU_starts a vector delineating starting indices of different reporting units in RU_vec
#' @param RU_vec a vector of collection indices, organized by reporting unit
#'
#' @examples
#' params <- tcf2param_list(alewife, 15)
#' SL <- geno_logL(params) %>% exp() %>% apply(2, function(x) x/sum(x))
#' correct <- avg_coll2correctRU(SL, params$coll, params$RU_starts, params$RU_vec)
#'
#' @export
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

#' generate a matrix of the average (by reporting unit) scaled likelihood of fish assigned to a reporting unit
#'
#' @param SL a scaled likelihood matrix
#' @param coll a vector of the collection of origin indices of the fish (length == ncol(SL))
#' @param RU_starts a vector delineating starting indices of different reporting units in RU_vec
#' @param RU_vec a vector of collection indices, organized by reporting unit
#'
#' @examples
#' params <- tcf2param_list(alewife, 15)
#' SL <- geno_logL(params) %>% exp() %>% apply(2, function(x) x/sum(x))
#' avg_SL_ru <- avg_assign_ru(SL, params$coll, params$RU_starts, params$RU_vec)
#'
#' @export
avg_assign_ru <- function(SL, coll, RU_starts, RU_vec) {
  ru_SL <- lapply(1:(length(RU_starts)-1), function(ru){
    apply(SL[RU_vec[(RU_starts[ru]+1):RU_starts[ru + 1]], ], 2, sum)
    }) %>%
    do.call("rbind", .)
  avg_ru_SL <- lapply(1:(length(RU_starts) - 1), function(ru) {
    as.array(ru_SL[, (coll %in% RU_vec[(RU_starts[ru]+1):RU_starts[ru + 1]]) ])
  }) %>%
    lapply(function(x) apply(x,MARGIN = 1, FUN = mean)) %>%
    simplify2array()
}

#' generate a matrix of the average (by collection) scaled likelihood of fish assigned to a reporting unit
#'
#' @param SL a scaled likelihood matrix
#' @param coll a vector of the collection of origin indices of the fish (length == ncol(SL))
#' @param RU_starts a vector delineating starting indices of different reporting units in RU_vec
#' @param RU_vec a vector of collection indices, organized by reporting unit
#'
#' @examples
#' params <- tcf2param_list(alewife, 15)
#' SL <- geno_logL(params) %>% exp() %>% apply(2, function(x) x/sum(x))
#' avg_SL_coll2ru <- avg_assign_coll(SL, params$coll, params$RU_starts, params$RU_vec)
#'
#' @export
avg_assign_coll <- function(SL, coll, RU_starts, RU_vec) {
  ru_SL <- lapply(1:(length(RU_starts)-1), function(ru){
    apply(SL[RU_vec[(RU_starts[ru]+1):RU_starts[ru + 1]], ], 2, sum)
  }) %>%
    do.call("rbind", .)
  avg_SL <- lapply(1:nrow(SL), function(x) as.array(ru_SL[, (coll == x)])) %>%
    lapply(function(x) apply(x,MARGIN = 1, FUN = mean)) %>%
    simplify2array()
}

