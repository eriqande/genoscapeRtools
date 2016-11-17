

#' read a set of 012 files from vcftools into a matrix of genotypes with appropriate dimnames
#'
#' Uses data.table::fread for speed. Note that the individual names should be unique! So should the
#' locus names.  The .pos file should have two columns of names (at least things are hard-wired for
#' that at the moment).
#' @param prefix  the (path) prefix on the files to be read.  The files will be \code{prefix.012},
#' \code{prefix.012.indv}, and \code{prefix.012.pos}
#' @param gz if this is set to TRUE, then the function expects that the 012 file ends in 012.gz
#' and is gzipped.  For this to work, you have to have a decent shell on your system with \code{gzcat}
#' for gunzipping the file to a stream.  If you are on a Mac, that should be fine.  If you are on
#' Windows, probably not.  If there are any problems, don't keep your file gzipped!
#' @return Returns a matrix with n-indiv rows and n-loci columns.  The rownames are the indv names and the
#' colnames are the pos's.
#' @export
read_012 <- function(prefix, gz = FALSE) {
  if(gz == TRUE) {
    file <- paste0("gzcat ", prefix, ".012.gz")
  } else {
    file <- paste0(prefix, ".012")
  }

  posfile <-  paste0(prefix, ".012.pos")
  indvfile <-  paste0(prefix, ".012.indv")

  ret <- data.table::fread(file, data.table = FALSE)[,-1] %>%
    as.matrix()

  posmat <- scan(posfile, what = "character") %>%
    matrix(ncol = 2, byrow = TRUE)
  pos <- paste(posmat[,1], posmat[,2], sep = "--")

  indv <- scan(indvfile, what = "character")

  dimnames(ret) <- list(indv = indv, pos = pos)
  ret
}
