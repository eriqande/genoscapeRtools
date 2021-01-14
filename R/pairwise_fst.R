
#' compute pairwise Fst using SNPRelate
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
pairwise_fst <- function(dat012, sample_groups) {
  new012 <- dat012 # just convenient in case we want to modify dat012 in a further update

  gds_file <- tempfile()

  SNPRelate::snpgdsCreateGeno(gds.fn = gds_file,
                              genmat = new012,
                              sample.id = rownames(new012),
                              snp.id = colnames(new012),
                              snpfirstdim = FALSE)

  gds_conn <- SNPRelate::snpgdsOpen(gds_file)  # open a connection to the gds file


  # to get all the pairwise combos of group in sample_groups
  # we use expand.grid
  sample_groups2 <- sample_groups %>%
    dplyr::filter(!is.na(group))

  pop_pairs <- expand.grid(pop1 = unique(sample_groups2$group),
                           pop2 = unique(sample_groups2$group),
                           stringsAsFactors = FALSE) %>%
    dplyr::filter(pop1 < pop2) %>%
    dplyr::arrange(pop1, pop2)


  # now we cycle over all those and get the Fst for each.
  # we do this as a summarise
  fst_calc_console_output <- capture.output(
    fst_data_frame <- pop_pairs %>%
      dplyr::group_by(pop1, pop2) %>%
      dplyr::summarise(Fst = computeFst(gds_conn, pop1, pop2, sample_groups2)) %>%
      dplyr::ungroup()
  )

  # then close the connection
  SNPRelate::snpgdsClose(gds_conn)

  # now, we clean up the output just a little
  # this is a hack.  The new version of SNPrelate puts this line into the output
  # that screws things up.  I am just going to remove it:
  fst_calc_console_output_stripped <- fst_calc_console_output[fst_calc_console_output != "Fst estimation on genotypes:"]

  console <- tibble::as_tibble(matrix(fst_calc_console_output_stripped, byrow = TRUE, ncol = 7))

  # and we grab the sample sizes out of it
  samp_sizes <- console %>%
    dplyr::mutate(tmp = stringr::str_replace_all(V6, "[,()]", "")) %>%
    dplyr::select(tmp) %>%
    dplyr::mutate(tmp = stringr::str_replace_all(tmp, "^ *", "")) %>%
    tidyr::separate(tmp, into = c("pop1", "n_pop1", "pop2", "n_pop2"), sep = " ") %>%
    dplyr::select(pop1, pop2, n_pop1, n_pop2)

  # and we can put the sample sizes on there
  fst_ret <- fst_data_frame %>%
    dplyr::left_join(samp_sizes, by = c("pop1", "pop2"))

  # then return a list
  list(
    Fst = fst_ret,
    stdout = console
  )
}



# a quick function to compute Fst given p1 and p2 an open  gds file and rn
# here rn is a tibble with a column sample and a column group
computeFst <- function(gds, p1, p2, rn) {
  vecs <- rn %>%
    dplyr::filter(group %in% c(p1, p2))

  ret <- SNPRelate::snpgdsFst(gds,
                              sample.id = vecs$sample,
                              population = as.factor(vecs$group),
                              method = "W&C84")

  ret$Fst
}
