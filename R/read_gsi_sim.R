

#' read a gsi_sim formatted input file into a tibble that rubias can use
#'
#' Note that this relies on a system call to awk.  It probably won't work
#' on Windows.
#' @param path path to the gsi_sim file
#' @param sample_type  should be "reference" or "mixture" depending on what kind of file it is
#' @param repunits the gsi_sim reporting units file.  Has no effect if sample_type is "mixture".  If
#' sample_type is "reference" and this is left as NULL, then all collections will be put in the the
#' "default_repu".
#' @export
read_gsi_sim <- function(path, sample_type, repunits = NULL) {

  tempd <- tempdir()
  dir.create(tempd, showWarnings = FALSE)
  gfile <- file.path(tempd, "gsifile.txt")
  rfile <- file.path(tempd, "repufile.txt")

  gscript <- system.file("scripts/clean-up-gsisim.awk", package = "rubias")
  rscript <- system.file("scripts/clean-up-gsisim-repunits.awk", package = "rubias")

  system(paste("awk -f", gscript, path, " > ", gfile))
  gsi_df <- readr::read_table2(gfile, na = "0")

  if (sample_type == "reference") {
    if (!is.null(repunits)) {
      system(paste("awk -f", rscript, repunits, " > ", rfile))
      repu_df <- readr::read_table2(rfile)
    } else {
      repu_df <- tibble::tibble(repunit = "default_repu", collection = unique(gsi_df$collection))
    }
    repu2 <- repu_df %>%
      dplyr::mutate(sample_type = "reference") %>%
      dplyr::select(sample_type, dplyr::everything())
  } else if (sample_type == "mixture") {
    repu_df <- tibble::tibble(repunit = NA_character_, collection = unique(gsi_df$collection))
    repu2 <- repu_df %>%
      dplyr::mutate(sample_type = "mixture") %>%
      dplyr::select(sample_type, dplyr::everything())
  } else {
    stop("sample_type must be either reference or mixture, not:  ", sample_type)
  }

  dplyr::left_join(repu2, gsi_df, by = "collection")

}
