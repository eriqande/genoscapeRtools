# functions related to investigating missing data in the 012 files and
# choosing suitable cutoff or index points for dropping/retaining individuals
# and loci


#' Order a genotype matrix by individual missingness at different locus numbers and plot results
#'
#' By default this just orders things and returns a data frame of missingness and a ggplot that can
#' be printed.  If you give it either \code{clean_pos} or \code{clean_indv} (or both) then an 012 file and indv
#' and pos files
#' will be written to  files named like "clean_file_prefix_indvXX_posYY" where XX is clean_indv and YY is clean_pos.
#' The order of pos and inv in that file will be the same as in the x012 input (except for the pos and indv that have
#' been removed, of course.) WARNING: If you positions have "--" in them, you will get some wonky results.
#' @param x012 a matrix of 012 genotypes with missing data coded as -1
#' @param reorder Logical.  if TRUE then the order of individuals is re-done for every number of loci
#' @param loc_nums vector of numbers of loci to investigate.  If not given, then it just does 10 equally-spaced
#' values.
#' @param clean_pos the number of loci to retain in a data set that has been "cleaned up".  Must be a single integer. If this is given, then
#' a cleaned set of 012 files will be written out keeping clean_pos loci.  If clean_indv is given, but this isn't, then all positions are
#' retained in the file output.
#' @param clean_indv the number of individuals to retain in a data set that has been "cleaned up". Must be a single integer. If this is not given,
#' but clean_pos is, then all the indv's will be returned.
#' @param clean_file_prefix the prefix to add to the output file if you either of retainL or retainI is given.
#' @param indv_specific_clean_pos Flag that I will probably eventually toss. I just have it in there for illustration.  If this is TRUE, when you are choosing
#' clean_pos markers, it will find the best ones for the given clean_indv you are retaining.  If FALSE, then it just takes the loci as ordered by missingness
#' amongst all individauls, which can lead to some loci having more missing data than you want!
#' @export
miss_curves_indv <- function(x012, reorder = TRUE, loc_nums = integer(0), clean_pos = integer(0), clean_indv = integer(0), clean_file_prefix = "cleaned", indv_specific_clean_pos = TRUE) {

  # default
  do_clean <- FALSE
  retain_indv_df <- NULL

  # deal with numbers of individuals and whether we will retain and clean some
  stopifnot(length(clean_indv) <= 1)
  stopifnot(length(clean_pos) <= 1)

  if(length(loc_nums) == 0) {
    loc_nums = floor(seq(1, ncol(x012), length.out = 11)[-1])
  }
  if(length(clean_indv) == 1) {
    do_clean <- TRUE
  }
  if(length(clean_pos) == 1) {
    do_clean <- TRUE
    loc_nums <- unique(sort(c(clean_pos, loc_nums)))  # add the clean_pos to the loc_nums
  }
  if(do_clean == TRUE) {
    if(length(clean_pos) == 0) {
      clean_pos <- ncol(x012)
    }
    if(length(clean_indv) == 0) {
      clean_indv <- nrow(x012)
    }
  }

  xmiss <- x012 == -1

  # compute the cumulative num of missing values
  locmiss <- colSums(xmiss)
  imiss <- rowSums(xmiss)
  xmiss_ord <- xmiss[order(imiss), order(locmiss)]
  cumul_miss <- mat_cumul_R(xmiss_ord, 1)   # This calls RCPP to avoid a slow ==> t(apply(xmiss_ord, 1, cumsum))


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
  if(reorder == TRUE || do_clean == TRUE) {
    smallfract_df <- smallfract_df %>%
      dplyr::group_by(NumLociKept) %>%
      dplyr::arrange(dplyr::desc(fract_non_missing)) %>%
      dplyr::mutate(xpos = 1:n())
  }

  # now, if we are writing out and returning a cleaned set, let's do that
  if(do_clean == TRUE) {
    fprefix <- paste(clean_file_prefix, "_indv", clean_indv, "_pos", clean_pos, ".012", sep = "")

    message("Picking out clean_pos and clean_indv and writing them to ", fprefix)
    retain_indv_df <- smallfract_df %>%
      dplyr::ungroup() %>%
      dplyr::filter(NumLociKept == clean_pos, xpos <= clean_indv)
    retain_indv <- retain_indv_df %>% dplyr::select(indv) %>% unlist() %>% unname()

    # now, get the positions to keep.  For this, we will use the best loci given only the retained
    # individuals
    if(indv_specific_clean_pos == TRUE) {
      nMiss <- sort(colSums(xmiss[retain_indv, ]))
      retain_pos <- names(nMiss)[1:clean_pos]
    } else {
      retain_pos <- colnames(xmiss_ord)[1:clean_pos]   # this is the lazy way I used to do it, which meant there were some sites that might have more missing data amongst the retained individauls than desired/expected
    }

    # now we pick those out in the original order
    clean_mat <- x012[rownames(x012) %in% retain_indv, colnames(x012) %in% retain_pos]

    cat(rownames(clean_mat), sep = "\n", file = paste(fprefix, ".indv", sep = ""))
    posnames <- colnames(clean_mat)  # get position names.  Then turn "--" into a \t

    cat(stringr::str_replace(posnames, pattern = "--", "\t"), sep = "\n", file = paste(fprefix, ".pos", sep = ""))
    rownames(clean_mat) <- 0:(nrow(clean_mat) - 1)
    write.table(clean_mat, file = fprefix, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

  }

  # now, prepare a ggplot of that as well
  g <- ggplot2::ggplot(smallfract_df, ggplot2::aes(x = xpos, y = fract_non_missing, colour = NumLocKeptF)) +
    ggplot2::geom_line() +
    ggplot2::xlab("Number of individuals (ordered by missingness)") +
    ggplot2::ylab("Fraction of non-missing genotypes in the individual")

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
    ggplot2::geom_line() +
    ggplot2::xlab("Number of loci (ordered by missingness)") +
    ggplot2::ylab("Fraction of non-missing genotypes at the locus")

  list(df = smallfract_df, plot = g)
}
