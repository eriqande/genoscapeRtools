
#' from on 012 file compute expected and observed genotype counts by population
#'
#' @param g012 an 012 matrix.  Missing data can be -1 or NA.
#' @param pops a data frame with one column "pop" and the other "sample" that tells
#' which populations individuals are in.  If "pop" entry is NA for an individual,
#' that individual is dropped silently. If NULL then everyone is assumed to be in the
#' same population.
#' @export
geno_freq_calcs <- function(g012, pops = NULL) {
  if (is.null(pops)) {
    return(geno_freq_calc_single(g012))
  }

  # need to add checking for input errors, like samples not in g012, etc
  popf <- dplyr::filter(pops, !is.na(pop))

  pop_list <- split(popf$sample, popf$pop)

  lapply(pop_list, function(x) geno_freq_calc_single(g012[x, ])) %>%
    dplyr::bind_rows(.id = "pop")
}



#' from on 012 file compute expected and observed genotype counts
#'
#' ...assuming HW equilbrium
#' @param g012 an 012 matrix.  Missing data can be -1 or NA.
geno_freq_calc_single <- function(g012) {
  g012[g012 == -1] <- NA

  gc <- 2 * colSums(!is.na(g012), na.rm = TRUE)  # number of non-missing gene copies
  p <- colSums(g012, na.rm = TRUE) / gc
  q <- 1.0 - p

  n0 <- colSums(g012 == 0, na.rm = TRUE)
  n1 <- colSums(g012 == 1, na.rm = TRUE)
  n2 <- colSums(g012 == 2, na.rm = TRUE)
  n <- n0 + n1 + n2

  tibble::tibble(snp = colnames(g012),
                 p = p,
                 ntot = n,
                 `0` = n0,
                 `1` = n1,
                 `2` = n2
                 ) %>%
    tidyr::gather(key = "geno", value = "n_obs", `0`, `1`, `2`) %>%
    dplyr::arrange(snp, geno) %>%
    dplyr::mutate(p_exp = ifelse(geno == 0, (1 - p)^2,
                               ifelse(geno == 1, 2 * p * (1 - p),
                                      ifelse(geno == 2, p^2, NA))),
                  n_exp = ntot * p_exp,
                  p_obs = n_obs / ntot,
                  z_score = (n_obs - n_exp) / sqrt(ntot * p_exp * (1 - p_exp))) %>%
    dplyr::select(snp, p, ntot, geno, p_exp, p_obs, n_exp, n_obs, z_score)


}
