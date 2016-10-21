

#' read a set of 012 files from vcftools into a matrix of genotypes with appropriate dimnames
#'
#' Uses data.table::fread for speed. Note that the individual names should be unique! So should the
#' locus names.  The .pos file should have two columns of names (at least things are hard-wired for
#' that at the moment).
#' @param prefix  the (path) prefix on the files to be read.  The files will be \code{prefix.012},
#' \code{prefix.012.indv}, and \code{prefix.012.pos}
#' @return Returns a matrix with n-indiv rows and n-loci columns.  The rownames are the indv names and the
#' colnames are the pos's.
#' @export
read_012 <- function(prefix) {
  file <- paste0(prefix, ".012")
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
