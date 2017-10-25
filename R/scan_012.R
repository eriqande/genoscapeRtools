#' read a set of 012 files from vcftools into a matrix of genotypes with appropriate dimnames
#'
#' This one skips data.table::fread and just scans things directly.
#' Note that the individual names should be unique! So should the
#' locus names.  The .pos file should have two columns of names (at least things are hard-wired for
#' that at the moment).
#' @param prefix  the (path) prefix on the files to be read.  The files will be \code{prefix.012},
#' \code{prefix.012.indv}, and \code{prefix.012.pos}
#' @param gz if this is set to TRUE, then the function expects that the 012 file ends in 012.gz
#' and is gzipped.  For this to work, you have to have a decent shell on your system with \code{gzcat}
#' for gunzipping the file to a stream.  If you are on a Mac, that should be fine.  If you are on
#' Windows, probably not.  If there are any problems, don't keep your file gzipped!
#' @param gzpos  logical flag saying whether the pos file is gzipped.
#' @param posfile_columns The number of columns in the .pos file.  This seems to usually be 1 or 2.
#' If this is 1, then it is assumed that the locus names are unique.  If it is greater than 1, then
#' the first two columns are catenated together with "--" to form the position names.
#' @return Returns a matrix with n-indiv rows and n-loci columns.  The rownames are the indv names and the
#' colnames are the pos's.
#' @export
scan_012 <- function(prefix, gz = FALSE, gzpos = FALSE, posfile_columns) {

  stopifnot(posfile_columns > 0)

  if (gz == TRUE) {
    file <- paste0(prefix, ".012.gz")
  } else {
    file <- paste0(prefix, ".012")
  }


  if (gzpos == TRUE) {
    posfile <-  paste0(prefix, ".012.pos.gz")
  } else {
    posfile <-  paste0(prefix, ".012.pos")
  }


  indvfile <-  paste0(prefix, ".012.indv")

  # first get the position names (and hence the number of them)
  message("Reading the position names ")
  pos_str <- matrix(scan(posfile, what = "character"),
                    ncol = posfile_columns,
                    byrow = TRUE)
  message("Read in ", length(pos_str), " positions.")
  if (posfile_columns == 1) {
    pos <- pos_str[,1]
  } else {
    pos <- paste0(pos_str[,1], "--", pos_str[,2])
  }
  message("First 10 position names are: ", paste(pos[1:10], collapse = ", "))


  # now, get the individual names
  message("Reading the individual names ")
  indv <- scan(indvfile, what = "character")


  message("Reading 012 matrix ")
  ret <- matrix(scan(file, what = integer()),
                byrow = TRUE,
                ncol = length(pos) + 1)[,-1]   # there is a extra column in the 012 file giving the index of the indv


  dimnames(ret) <- list(indv = indv, pos = pos)
  ret
}
