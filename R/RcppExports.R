# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Compute average pairwise nucleotide diversity statistics for positions from a to b
#'
#' This is the internal Rcpp function that will be called once for each window.
#' @param x an 0,1,2,-1 integer matrix
#' @param a the lowest (base-1) index of the SNPs in the window
#' @param b the highest (base-1) index of the SNPs in the window
#' @return This returns a list of vectors.  First vector is the average pairwise nucleotide difference
#' for each pair of individuals.  Second vector is the fraction of sites missing in the pair.
#' @export
apwnd_window_internal <- function(x, a, b) {
    .Call('_genoscapeRtools_apwnd_window_internal', PACKAGE = 'genoscapeRtools', x, a, b)
}

#' Sample a genotype (0,1,2) from probabilities for all the snps and individuals
#' in x.
#' @param x a numeric vector that holds the posterior probs from an
#' ANGSD geno.gz file (from doGenos 32)
#' @param samples character vector of sample names
#' @param snps character vector of SNP names
#' @param thresh any snp in an individual must have one posterior greater than T in order to
#' have a non-missing value simulated for it.  For example, if T is 0.6 and the posteriors
#' are 0.33, 0.33, 0.33, this locus will just get a missing (-1).
#' @export
sample_from_angsd_probs <- function(x, samples, snps, thresh) {
    .Call('_genoscapeRtools_sample_from_angsd_probs', PACKAGE = 'genoscapeRtools', x, samples, snps, thresh)
}

#' return a matrix of cumumlative sums along the rows or colums
#' @param x a numeric matrix
#' @param dim 1 for rows, 2 for columns
mat_cumul_cpp <- function(x, dim) {
    .Call('_genoscapeRtools_mat_cumul_cpp', PACKAGE = 'genoscapeRtools', x, dim)
}

