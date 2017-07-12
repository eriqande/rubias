

#' Write a mixerence data frame to gsi_sim format baseline and repunits file
#'
#' Note, this is only intended to work with integer-valued alleles, at the moment.
#' It was just written for testing and verifying that things are working correctly.
#'
#' @param mix mixture data frame
#' @param gen_start_col column in which the genetic data start
#' @param mixprefix path to write the mixture file to.  The mixture collection name + .txt will
#' be appended to this.
#' @export
write_gsi_sim_mixture <- function(mix, gen_start_col, mixprefix = "mixture") {

  locus_names <- names(mix)[seq(gen_start_col, ncol(mix), by = 2)]
  loccols <- names(mix)[gen_start_col:ncol(mix)]

  # then write the full mixture
  mix[is.na(mix)] <- 0
  mix_list <- split(mix, mix$collection)

  dump <- lapply(names(mix_list), function(x) {

    mixout <- paste(mixprefix, "-", x, ".txt", sep = "")
    cat(nrow(mix), (ncol(mix) - gen_start_col + 1) / 2, "\n", file = mixout)  # number of indivs and loci on top line
    cat(locus_names, sep = "\n", file = mixout, append = TRUE)
    cat("POP", x, "\n", file = mixout, append = TRUE)
    write.table(mix_list[[x]][, c("indiv", loccols)], sep = "  ", file = mixout, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  })

}
