% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hwe_exact_test.R
\name{testHWE}
\alias{testHWE}
\title{a quick function to do the HW exact test on on an open  gds file}
\usage{
testHWE(gds, sample_ids__, d012)
}
\arguments{
\item{gds}{the gds connection}

\item{sample_ids__}{a vector of sample IDs}

\item{d012}{the original 012 file}
}
\description{
Uses the individual ids in sample_ids__
this returns a tibble so that it can be used with dplyr::do().
It also returns the observed genotype calcs for each locus.
}
