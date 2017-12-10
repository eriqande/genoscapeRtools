
#' compute fraction of heterozygous and homozygous and missing SNPs in a sliding window
#'
#' This is just a quick thing Eric threw together.  It takes an 012,-1 matrix in which
#' the rownames are the sample IDs and the column names are the positions (without any
#' chromosome name).  Then it does kernel smoothing to estimate the fraction of SNPs
#' in each chunk that are homozygous reference, homozygous alternate, heterozygous, and
#' missing. These all get returned in a tidy tibble.
#'
#' This "_1" version is intended to operate on a single chromosome.
#'
#' @param d012 the 012,-1 matrix with column names being the positions
#' @param width the sliding window width.  Each window will be a non-overlapping
#' window of this width, and the position of it will be reported as the center of
#' the window. Default width is 10 Kb.
#' @export
#' @example
#' dat <- scan_012("~/Documents/git-repos/omy28-investigations/wgs_play/more_big_region_wgs", gz = FALSE, posfile_columns = 2)
#' colnames(dat) <- stringr::str_replace(colnames(dat), "omy28--", "")
hohz_windows_1 <- function(d012, width = 1e04) {

  # positions of the SNPs as numeric
  pts <- as.numeric(colnames(d012))

  # center points of the windows
  x.points <- seq(from = pts[1] + 0.5 * width,
                  to = pts[length(pts)],
                  by = width)

  zero_mat <- apply(d012 == 0, 1, function(y) ksmooth(x = pts, y = y, kernel = "box", bandwidth = width, x.points = x.points)$y)
  one_mat <- apply(d012 == 1, 1, function(y) ksmooth(x = pts, y = y, kernel = "box", bandwidth = width, x.points = x.points)$y)
  two_mat <- apply(d012 == 2, 1, function(y) ksmooth(x = pts, y = y, kernel = "box", bandwidth = width, x.points = x.points)$y)
  miss_mat <- apply(d012 == -1, 1, function(y) ksmooth(x = pts, y = y, kernel = "box", bandwidth = width, x.points = x.points)$y)


  # now, we also want to count up the number of SNPs per base-pair in each of those windows
  bin_endpts <- c(pts[1] - 0.001, x.points + 0.5 * width)  # bins corresponding to the windows
  bin_endpts[length(bin_endpts)] <- pts[length(pts)] + 0.001  # the last one has to be corrected
  bin_widths <- bin_endpts[-1] - bin_endpts[-length(bin_endpts)]  # the widths of each of those bins
  snps_per_bp <- as.vector(table(cut(pts, breaks = bin_endpts))) / bin_widths

  # now return these as tidy things
  ret1 <- tibble::tibble(sample = rep(rownames(d012), each = nrow(zero_mat)),
                         pos = rep(x.points, times = ncol(zero_mat)),
                         hom_ref = as.vector(zero_mat),
                         het = as.vector(one_mat),
                         hom_alt = as.vector(two_mat),
                         missing = as.vector(miss_mat),
                         snps_per_bp = rep(snps_per_bp, times = ncol(zero_mat)))

  ret1

}
