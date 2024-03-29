% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/explore-missingness.R
\name{miss_curves_indv}
\alias{miss_curves_indv}
\title{Order a genotype matrix by individual missingness at different locus numbers and plot results}
\usage{
miss_curves_indv(
  x012,
  reorder = TRUE,
  loc_nums = integer(0),
  clean_pos = integer(0),
  clean_indv = integer(0),
  clean_file_prefix = "cleaned"
)
}
\arguments{
\item{x012}{a matrix of 012 genotypes with missing data coded as -1}

\item{reorder}{Logical.  if TRUE then the order of individuals is re-done for every number of loci}

\item{loc_nums}{vector of numbers of loci to investigate.  If not given, then it just does 10 equally-spaced
values.}

\item{clean_pos}{the number of loci to retain in a data set that has been "cleaned up".  Must be a single integer. If this is given, then
a cleaned set of 012 files will be written out keeping clean_pos loci.  If clean_indv is given, but this isn't, then all positions are
retained in the file output.}

\item{clean_indv}{the number of individuals to retain in a data set that has been "cleaned up". Must be a single integer. If this is not given,
but clean_pos is, then all the indv's will be returned.}

\item{clean_file_prefix}{the prefix to add to the output file if you either of retainL or retainI is given.}
}
\description{
By default this just orders things and returns a data frame of missingness and a ggplot that can
be printed.  If you give it either \code{clean_pos} or \code{clean_indv} (or both) then an 012 file and indv
and pos files
will be written to  files named like "clean_file_prefix_indvXX_posYY" where XX is clean_indv and YY is clean_pos.
The order of pos and inv in that file will be the same as in the x012 input (except for the pos and indv that have
been removed, of course.) WARNING: If you positions have "--" in them, you will get some wonky results.
}
