# functions related to investigating missing data in the 012 files and
# choosing suitable cutoff or index points for dropping/retaining individuals
# and loci


#' Order a genotype matrix by individual missingness at different locus numbers and plot results
#'
#' @param x012 a matrix of 012 genotypes with missing data coded as -1
#' @param reorder Logical.  if TRUE then the order of individuals is re-done for every number of loci
#' @param loc_nums vector of numbers of loci to investigate.  If not given, then it just does 10 equally-spaced
#' values.
#' @export
miss_curves_indv <- function(x012, reorder = TRUE, loc_nums = integer(0)) {
  xmiss <- x012 == -1

  # compute the cumulative num of missing values
  locmiss <- colSums(xmiss)
  imiss <- rowSums(xmiss)
  xmiss_ord <- xmiss[order(imiss), order(locmiss)]
  cumul_miss <- mat_cumul_R(xmiss_ord, 1)   # This calls RCPP to avoid a slow ==> t(apply(xmiss_ord, 1, cumsum))

  if(length(loc_nums) == 0) {
    loc_nums = floor(seq(1, ncol(x012), length.out = 11)[-1])
  }

  # select those loci and then make a data frame of it
  small <- cumul_miss[, loc_nums]
  smallfract <- small / rep(loc_nums, each = nrow(small))
  smallfract_df <- as.data.frame(smallfract) %>%
    dplyr::tbl_df() %>%
    setNames(loc_nums) %>%
    dplyr::mutate(xpos = 1:n(),
                  indv = rownames(small)) %>%
    tidyr::gather(data = ., key = "NumLociKept", value = "fract_missing", convert = TRUE, -xpos, -indv) %>%
    dplyr::mutate(NumLocKeptF = factor(NumLociKept),
           fract_non_missing = 1 - fract_missing)


  # now, if we are to reorder it, we can just modify the xpos column grouped by NumLociKept
  if(reorder == TRUE) {
    smallfract_df <- smallfract_df %>%
      dplyr::group_by(NumLociKept) %>%
      dplyr::arrange(dplyr::desc(fract_non_missing)) %>%
      dplyr::mutate(xpos = 1:n())
  }

  # now, prepare a ggplot of that as well
  g <- ggplot2::ggplot(smallfract_df, ggplot2::aes(x = xpos, y = fract_non_missing, colour = NumLocKeptF)) +
    ggplot2::geom_line()

  list(df = smallfract_df, plot = g)
}





#' Order a genotype matrix by locus missingness at different numbers of individuals and plot results
#'
#' @param x012 a matrix of 012 genotypes with missing data coded as -1
#' @param reorder Logical.  if TRUE then the order of loci is re-done for every number of loci
#' @param locus_thin_to For speed with plotting it is possible to retain just  fraction of the loci,
#' once the cumulatives have been counted.  Default is 1000 of them.
#' @param indv_nums The number of individuals to explore retaining.  By default it just does them
#' in 10 roughly equal chunks.
#' @export
miss_curves_locus <- function(x012, reorder = TRUE, locus_thin_to = 5000, indv_nums = integer(0)) {
  xmiss <- x012 == -1

  # compute the cumulative sum of missing values
  locmiss <- colSums(xmiss)
  imiss <- rowSums(xmiss)
  xmiss_ord <- xmiss[order(imiss), order(locmiss)]
  cumul_miss <- mat_cumul_R(xmiss_ord, 2)


  if(length(indv_nums) == 0) {
    indv_nums = floor(seq(1, nrow(x012), length.out = 11)[-1])
  }
  thin <- floor(seq(1,ncol(x012), length.out = locus_thin_to))

  # select those loci and then make a data frame of it
  small <- cumul_miss[indv_nums, thin]
  colnames(small) <- thin
  rownames(small) <- indv_nums
  smallfract <- small / rep(indv_nums, nrow(small))
  smallfract_df <- as.data.frame(smallfract) %>%
    dplyr::tbl_df() %>%
    dplyr::mutate(num_indv = factor(indv_nums)) %>%
    tidyr::gather(data = ., key = "num_loci", value = "fract_missing", convert = TRUE, -num_indv) %>%
    dplyr::mutate(fract_non_missing = 1 - fract_missing) %>%
    dplyr::group_by(num_indv) %>%
    dplyr::arrange(dplyr::desc(fract_non_missing)) %>%
    dplyr::mutate(num_loci_reordered = sort(num_loci)) %>%
    dplyr::ungroup()




  # now, prepare a ggplot of that as well
  g <- ggplot2::ggplot(smallfract_df, ggplot2::aes(x = num_loci_reordered, y = fract_non_missing, colour = num_indv)) +
    ggplot2::geom_line()

  list(df = smallfract_df, plot = g)
}
