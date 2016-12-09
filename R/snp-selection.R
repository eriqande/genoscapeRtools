



#' return SNPs and say whether the are assayable or not, along with the MAFs
#'
#' This should be run on a set of data that has not been too highly filtered---you
#' want to know if there is much variation, even if it wasn't seen in everyone.
#' @param w a matrix like produced by read_012
#' @param MAF the minimum minor allele frequency to retain a SNP either for assaying
#' or for determining whether there is flanking variation.
#' @param flank  If a SNP is to be assayable, it needs to have no SNPs within flank
#' base pairs of itself.
#' @return This returns a data frame with the columns:
#' \enumerate{
#' \item snp num_gene_copies
#' \item freq
#' \item scaffold
#' \item bp
#' \item length
#' \item leftpos
#' \item rightpos
#' \item leftpad
#' \item rightpad
#' \item assayable
#' \item on_right_edge
#' \item on_left_edge
#' }
#' @export
compile_assayable <- function(w, minMAF = 0.03, flank = 20) {
  w[w == -1] <- NA  # mark the missing ones as missing.
  gc <- (2 * colSums(!is.na(w)))  # count how many gene copies are copied at each SNP
  frq <- colSums(w, na.rm = TRUE) / gc  # get the  overall freqs of the "1" allele
  keepers <- frq > minMAF & frq < 1.0 - minMAF

  # now make a data frame of them and keep track of the scaffold and base pairs too:
  df <- dplyr::data_frame(snp = colnames(w)[keepers],
                          num_gene_copies = gc[keepers],
                          freq = frq[keepers]) %>%
    dplyr::mutate(pos = snp) %>%
    tidyr::separate(pos, into = c("scaffold", "bp"), sep = "--", convert = TRUE) %>%
    dplyr::mutate(length = as.numeric(str_replace(scaffold, "scaff.*size", ""))) %>%
    dplyr::group_by(scaffold) %>%
    dplyr::mutate(leftpos = c(0,bp)[-(n() + 1)],
                  rightpos = c(bp, length[1])[-1]) %>%
    dplyr::mutate(leftpad = bp - leftpos,
                  rightpad = rightpos - bp) %>%
    dplyr::ungroup()

  # and in the end, we will specify the status of each:
  df %>%
    dplyr::mutate(assayable = ifelse(rightpad > 20 & leftpad > 20, TRUE, FALSE)) %>%
    dplyr::mutate(on_right_edge = ifelse(rightpad > 2000, TRUE, FALSE),
                  on_left_edge = ifelse(leftpad > 2000, TRUE, FALSE)
    )
}




#' count up number of observed alleles in each group
#'
#' boing.
#' @param x012 an 012 matrix like that returned by read_012
#' @param glist a named list of ID vectors. Each name is the name of the group
#' and each component is a vector of ID names that correspond to some of the individuals
#' in the rownames of x012.
#' @return Returns a data frame with columns pos, n0, n1, group, ntot, and freq
#' @export
group_012_cnts <- function(x012, glist) {
  x012[x012 == -1] <- NA  # turn 1 to NA
  lapply(glist, function(x) {
    y <- x012[x,]
    ntot <- 2 * colSums(!is.na(y))
    n1 <- colSums(y, na.rm = TRUE)
    n0 <- ntot - n1
    tibble(pos = names(ntot), n0 = n0, n1 = n1)
  }) %>%
    bind_rows(.id = "group") %>%
    mutate(ntot = n0 + n1,
           freq = n1 / ntot)
}





#' from vectors of snp allele freqs, compute the expected prob of correct assignment
#'
#' this is most useful for use with dplyr and mutate
#' @param p0 allele frequency in first population
#' @param p1 allele frequency in second population
#' @return a vector of the same length as p0 and p1
#' @export
correct_ass_prob <- function(p0, p1) {
  # genotype freqs in the first population
  P00 <- (1 - p0)^2
  P12 <- 2 * p0 * (1 - p0)
  P11 <- p0^2

  # same in the second population
  Q00 <- (1 - p1)^2
  Q12 <- 2 * p1 * (1 - p1)
  Q11 <- p1^2

  # prob of correctly assigning a P-population individual:
  PCP <- P00 * (P00 >= Q00)    +   P12 * (P12 >= Q12)    +    P11 * (P11 >= Q11)

  # same for a Q-population individual
  PCQ <- Q00 * (P00 < Q00)    +   Q12 * (P12 < Q12)    +    Q11 * (P11 < Q11)

  # return the arithmetic mean of those two
  (PCP + PCQ) / 2

}




#' make a comparison between two groups that are in df
#'
#' boing
#' @param df the data frame you are picking things from.  This should be
#' like the data frame that you get back from group_012_cnts
#' @param grp1 the name of the first group in the comparison
#' @param grp2 the name of the second group in the comparison
#' @return Sneds back a data frame sorted in descending order of
#' correct assignment prob.
#' @export
comp_pair <- function(df, grp1, grp2) {
  x <- df %>%
    filter(group == grp1) %>%
    select(group, pos, freq, ntot)

  y <- df %>%
    filter(group == grp2) %>%
    select(group, pos, freq, ntot)

  inner_join(x, y, by = c("pos")) %>%
    mutate(comparison = paste(group.x, group.y, sep = " vs ")) %>%
    select(pos, comparison, starts_with("group"), starts_with("freq"), starts_with("ntot")) %>%
    mutate(cap = correct_ass_prob(freq.x, freq.y)) %>%
    arrange(desc(cap))
}
