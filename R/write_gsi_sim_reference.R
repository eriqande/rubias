

#' Write a reference data frame to gsi_sim format baseline and repunits file
#'
#' Note, this is only intended to work with integer-valued alleles, at the moment.
#' It was just written for testing and verifying that things are working correctly.
#'
#' @param  ref reference data frame
#' @param gen_start_col column in which the genetic data start
#' @param baseout path to write the baseline file to. Required.
#' @param repout path to write the repunits file to. Required.
#' @export
#' @examples
#' # create a temp directory to put example outputs
#' dd <- tempdir()
#' basefile <- file.path(dd, "baseline.txt")
#' repunitsfile <- file.path(dd, "repunits.txt")
#'
#' # print those
#' basefile
#' repunitsfile
#'
#' # note that in practice you will probably want to specify
#' # your own filepaths...
#'
#' # run the function
#' write_gsi_sim_reference(alewife, 17, basefile, repunitsfile)
write_gsi_sim_reference <- function(ref, gen_start_col, baseout = "baseline.txt", repout = "repunits.txt") {

  # first, write the reporting unit file
  reps_list <- ref %>%
    dplyr::count(repunit, collection) %>%
    dplyr::select(-n) %>%
    dplyr::arrange(repunit, collection) %>%
    base::split(.$repunit)

  if (file.exists(repout)) file.remove(repout)
  dump <- lapply(names(reps_list), function(x){
    cat("REPUNIT", x, "\n", file = repout, append = TRUE)
    cat(paste("    ", reps_list[[x]]$collection, "\n", sep = ""), sep = "", file = repout, append = TRUE)
  })


  # then write the full baseline
  ref[is.na(ref)] <- 0
  ref_list <- split(ref, ref$collection)

  cat(nrow(ref), (ncol(ref) - gen_start_col + 1) / 2, "\n", file = baseout)  # number of indivs and loci on top line
  locus_names <- names(ref)[seq(gen_start_col, ncol(ref), by = 2)]
  cat(locus_names, sep = "\n", file = baseout, append = TRUE)

  loccols <- names(ref)[gen_start_col:ncol(ref)]
  dump <- lapply(names(ref_list), function(x) {
    cat("POP", x, "\n", file = baseout, append = TRUE)
    write.table(ref_list[[x]][, c("indiv", loccols)], sep = "  ", file = baseout, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  })

}
