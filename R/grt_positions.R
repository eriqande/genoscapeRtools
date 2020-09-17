

#' extract the positions from the locus names in an 012 file
#'
#' Simple little function.  You pass it an 012 matrix that was read
#' in by \code{\link{scan_012}} for example, and it gives the positions
#' as an integer vector.  By default it assumes that the markers are named
#' "chrom--pos", but you can supply whatever regex will pick that out.
#' I think that I should have figured it out so that it this will work even
#' if the chrom name as a -- in it.
#' @param dat012 and 012 matrix
#' @param rgx the regex expression that will return the number part of the marker name.
#' This is set up to return the first capture group of the regex
#' @export
grt_positions <- function(dat012, rgx = "--([0-9]+)$") {
  as.integer(stringr::str_match(colnames(dat012), rgx)[,2])
}
