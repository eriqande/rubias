

#' check a baseline and mixture file together to ensure the known_collections are valid if they exist
#'
#' Simple function that checks for known_collections columns in a reference and mixture and
#' makes sure that they are compliant.
#' If there is a non-NA entry in the Mixture frame's known_collection column this
#' function returns TRUE.  Otherwise it returns FALSE.
#' @param R reference data frame
#' @param M mixture data frame
#' @export
#' @keywords internal
check_known_collections <- function(R, M) {

  ret <- FALSE

  if(!is.null(R[["known_collection"]]) && is.null(M[["known_collection"]])) {
    stop("Reference data frame has known_collection column, but not the mixture data frame.")
  }
  if(!is.null(M[["known_collection"]]) && is.null(R[["known_collection"]])) {
    stop("Mixture data frame has known_collection column, but not the reference data frame.")
  }

  if(!is.null(R[["known_collection"]])) {
    if(!is.character(R[["known_collection"]])) {
      stop("Column known_collection in reference data frame must be a character vector.")
    }
  }
  if(!is.null(M[["known_collection"]])) {
    if(!is.character(M[["known_collection"]])) {
      stop("Column known_collection in mixture data frame must be a character vector.")
    }
    KC <- M[["known_collection"]]
    uKC <- unique(KC[!is.na(KC)])
    mKC <- setdiff(uKC, unique(R[["collection"]]))
    if(length(mKC) > 0) {
      stop("These known_collection entries in mixture data frame not found amongst reference collections: ",
           paste(mKC, collapse = ", "), ". Please fix that and try again.")
    }
    ret <- any(!is.na(KC))
  }
  ret
}
