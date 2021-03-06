% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hohz_windows.R
\name{hohz_windows_1}
\alias{hohz_windows_1}
\title{compute fraction of heterozygous and homozygous and missing SNPs in a sliding window}
\usage{
hohz_windows_1(d012, width = 10000)
}
\arguments{
\item{d012}{the 012,-1 matrix with column names being the positions}

\item{width}{the sliding window width.  Each window will be a non-overlapping
window of this width, and the position of it will be reported as the center of
the window. Default width is 10 Kb.}
}
\description{
This is just a quick thing Eric threw together.  It takes an 012,-1 matrix in which
the rownames are the sample IDs and the column names are the positions (without any
chromosome name).  Then it does kernel smoothing to estimate the fraction of SNPs
in each chunk that are homozygous reference, homozygous alternate, heterozygous, and
missing. These all get returned in a tidy tibble.
}
\details{
This "_1" version is intended to operate on a single chromosome.
}
\examples{
dat <- scan_012("~/Documents/git-repos/omy28-investigations/wgs_play/more_big_region_wgs", gz = FALSE, posfile_columns = 2)
colnames(dat) <- stringr::str_replace(colnames(dat), "omy28--", "")
}
