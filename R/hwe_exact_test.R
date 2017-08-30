

#' perform exact tests of H-W disequilibrium for each locus using SNPRelate
#'
#' This just wraps SNPRelate's function and makes it easy to compute
#' Fst for all pairwise combinations.  Returns a tidy data frame.
#' @param dat012  An 012 matrix.  Missing data should be -1.
#' @param sample_groups this should be a tibble with a column "sample"
#' and another column "group".  The Fst calculations are done pairwise between
#' the groups.  You probably don't want individuals in sample_groups that are
#' not in the rownames of dat012.  If group is NA for any individuals, they
#' get silently tossed out.
#' @export
hwe_exact_test <- function(dat012, sample_groups) {
  new012 <- dat012 # just convenient in case we want to modify dat012 in a further update

  gds_file <- tempfile()

  SNPRelate::snpgdsCreateGeno(gds.fn = gds_file,
                              genmat = new012,
                              sample.id = rownames(new012),
                              snp.id = colnames(new012),
                              snpfirstdim = FALSE)

  gds_conn <- SNPRelate::snpgdsOpen(gds_file)  # open a connection to the gds file


  # cycle over the groups and return a tibble for each
  sample_groups2 <- sample_groups %>%
    dplyr::filter(!is.na(group))

  ret <- sample_groups2 %>%
    dplyr::group_by(group) %>%
    dplyr::do(testHWE(gds_conn, .$sample, new012)) %>%
    ungroup()


  # then close the connection
  SNPRelate::snpgdsClose(gds_conn)

  ret
}



#' a quick function to do the HW exact test on on an open  gds file
#'
#' Uses the individual ids in sample_ids__
#' this returns a tibble so that it can be used with dplyr::do().
#' It also returns the observed genotype calcs for each locus.
#' @param gds the gds connection
#' @param sample_ids__ a vector of sample IDs
#' @param d012 the original 012 file
testHWE <- function(gds, sample_ids__, d012) {


  ret <- SNPRelate::snpgdsHWE(gds,
                              sample.id = sample_ids__,
                              with.id = TRUE
  )

  pvals <- tibble::tibble(snp = ret$snp.id, p_value = ret$pvalue)

  geno_nums <- tibble::tibble(
    n0 = colSums(d012 == 0, na.rm = TRUE),  # add na.rm==TRUE just in case someone has changed -1 to NA
    n1 = colSums(d012 == 1, na.rm = TRUE),
    n2 = colSums(d012 == 2, na.rm = TRUE)
  )

  dplyr::bind_cols(pvals, geno_nums)
}
