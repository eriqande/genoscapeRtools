

#' break markers up into non-overlapping windows and compute ave-pw-nuc-div
#'
#' more needed
#' @param d012 an 012 matrix.  Positions must be strictly increasing (so, no more than one
#' chromosome / scaffold)
#' @param width the desired window width in base pairs
#' @param positions can be used to pass the positions in if they aren't correctly parsed
#' in the 012 file
#' @export
apwnd_windows <- function(d012, width, positions = NULL) {
  if (is.null(positions)) {
    positions = grt_positions(d012)
  }

  # find the windows.  I'm going to be really lax about the very last one---just going to
  # make it of width "width" even if there is not much data in that.
  maxp <- positions[length(positions)]
  minp <- positions[1]
  extra <- 0
  if (((maxp - minp) %% width) > 0) {
    extra <- 1
  }
  bounds <- seq(minp, width * ((maxp %/% width) + extra), by = width)

  # cut the positions into bins
  bins <- cut(positions, bounds)

  # then find the midoints of each bin
  mids <- ((bounds[-length(bounds)] + bounds[-1]) / 2)

  # get a list of starting and ending positions in each window
  windpos <- split(1:length(positions), bins)

  # now cycle over those positions and get the pairwise nds.  Some windows might be empty
  # and we return 0s for that
  boing <- lapply(windpos, function(x) {
    if (length(x) > 0) {
      apwnd_window_internal(d012, a = min(x), b = max(x))
    } else {
      list(nd = 0, msum = 0)
    }
  })

  names(boing) <- mids

  ret <- dplyr::bind_rows(boing, .id = "window_midpoint") %>%
    dplyr::mutate(window_midpoint = as.numeric(window_midpoint))

  ret
}
