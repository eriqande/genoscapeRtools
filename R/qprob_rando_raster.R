

#' return an RGB raster of colors to depict breeding pops on a map
#'
#' This takes the raster brick that comes out of tess3Q_map_rasters() and
#' does some rescaling and random sampling to return a raster of colors that
#' I hope will do a reliable job of representing (in some way) predicted assignment accuracy
#' over space.  This blends the colors in the intermediate areas by just sampling
#' pixels of different colors according to the smoothed Q values we get out of
#' tess3Q_map_rasters().
#' @param TRB the raster brick that comes out of tess3Q_map_rasters.  Typically you
#' will want to name the layers in it before you toss it into this function.
#' @param cols a named vector of colors. The names must correspond to the names of TRB.
#' @param alpha Transparency value from 0 (totally invisible) to 1 (totally opaque).
#' @param alpha_scale a number between 0 and 1 which gets multiplied by the scaled (so the max
#' is one and the min is the same) max-regions prob values.  This will override the alpha parameter
#' without a warning.
#' @param scale_min_to This lets you scale the minimum values for setting alpha so that in each max area they are
#' scale_at_min.
#' @param abs_thresh the smoothed Q value above which you always give the
#' cell the color it belongs to.  This is the absolute version, which means that
#' once the smoothed Q values have been scaled in a focal region, any of them
#' greater than abs_thresh will effectively be treated as 1s.
#' @param rel_thresh Similar to abs_thresh but it does it on a relative basis.  For
#' example, if you choose .80, then everything between the max and 80% of the way to the
#' min value in the max focal region will be assigned 1s.
#' @param alpha_chop_max After all the scaling and squishing, if any values of alpha (on the 0..255 scale)
#' are greater than this, then they get set to alpha_chop_max.
#' @export
#' @examples
#' data(wifl_stack)
#' wifl_colors <-  c(
#'  INW = "#984ea3",
#'  EST = "#377eb8",
#'  PNW = "#4daf4a",
#'  SSW = "#ff7f00",
#'  SCC = "#ffff33",
#'  KER = "#e41a1c",
#'  WMT = "#00ffff")
#'
#'  rr <- qprob_rando_raster(TRB = wifl_stack, cols = wifl_colors, alpha_scale = 2.0, abs_thresh = 0.0, alpha_exp = 1.45)
#'
qprob_rando_raster <- function(TRB, cols, alpha = 0.7, alpha_scale = NULL, scale_min_to = 0, abs_thresh = NULL, rel_thresh = NULL, alpha_exp = 1, alpha_chop_max = 255) {

  if(!is.null(abs_thresh) && !is.null(rel_thresh) ) stop("Only one of abs_thresh or rel_thresh can be specified.")

  screwup <- setdiff(names(cols), names(TRB))
  if(length(screwup) > 0) stop("Mismatch between names(TRB) and names(cols). Offenders: ", paste(screwup, collapse = ", "))

  # Now, make a stack that has values for each group, only where
  # that group has maximal values.
  max_groups <- TRB
  maxes <- max(max_groups)
  keep_cells <- max_groups >= maxes
  max_groups[keep_cells <= 0] <- NA

  # then create a scaled version, where each value within
  # a group-max area is scaled so the max value is 1.0
  # but the lowest value remains what it is.
  scal <- max_groups
  a <- cellStats(max_groups, min, na.rm = TRUE) # the min values in each of those spots
  b <- cellStats(max_groups, max, na.rm = TRUE) # the max values in each of those spots
  for (i in seq_along(a)) {
    scal[[i]] <-  max_groups[[i]] / (max_groups[[i]] * (b[i] - 1) / (b[i] - a[i])  + (b[i] - a[i] * b[i]) / (b[i] - a[i]))
  }
  names(scal) <- names(TRB)

  # Now, we are going to use those values to rescale everything else.
  # So convert NAs to 0
  scal_z <- scal
  scal_z[is.na(scal_z)] <- 0


  # now, insert those scaled values back in with the other
  # values outside of their max areas

  # force values < 0 to be 0
  TRB_z <- TRB
  TRB_z[TRB_z < 0] <- 0

  ## making scal_z_thresh, in which all values > thresh were turned into 1s, and then
  ## using that in place of scal_z.
  scal_z_thresh <- scal_z
  if(!is.null(abs_thresh)) {
    abt <- rep(abs_thresh, length.out = nlayers(scal_z))  # just doing this so that we could, if we wanted to, specify different cutoffs for different focal regions
    for(i in 1:nlayers(scal_z_thresh)) {
      scal_z_thresh[[i]][scal_z_thresh[[i]] > abt[i]] <- 1
    }
  }
  if(!is.null(rel_thresh)) {
    rel <- rep(rel_thresh, length.out = nlayers(scal_z))  # just doing this so that we could, if we wanted to, specify different cutoffs for different focal regions
    mins <- cellStats(scal, min, na.rm = TRUE)
    threshes <- 1 - (1 - mins) * 0.8
    for(i in 1:nlayers(scal_z_thresh)) {
      scal_z_thresh[[i]][scal_z_thresh[[i]] > threshes[i]] <- 1
    }
  }


  # give each area the scaled scal_z_thresh values where they are maximal
  # and the unscaled values outside of that area
  TRB_rescaled <- (TRB_z * !keep_cells) + scal_z_thresh



  # here we sum up the values of the probabilities in that focal max regions
  # _from the non-max cluster_
  for(i in 1:nlayers(TRB_rescaled)) {
    tmpbrick <- TRB_rescaled * keep_cells[[i]]  # picks out just the focal region
    tmpbrick[[i]] <- 0  # removes values of the focal cluster from the sum

    non_focal_sums <- sum(tmpbrick)  # this is a raster with the sums of non-focal clusters in the focal region
    non_focal_sums[non_focal_sums == 0] <- 1  # we will be dividing these sums out for these regions, and want the other regions to be unaffected, so we set those to 1s




    # now we rescale all those non-focal clusters
    for (j in 1:nlayers(TRB_rescaled)) {
      if (i != j) {
        TRB_rescaled[[j]] <- (1 - scal_z_thresh[[i]]) * TRB_rescaled[[j]] / non_focal_sums
      }
    }
  }


  # at the end of that we have some tiny values that are negative.  Let's turn those to 0s
  TRB_rescaled[TRB_rescaled < 0] <- 0

  names(TRB_rescaled) <- names(TRB)

  # all right! Now, let's sample the cluster from these probabilities, so that
  # we get a raster in which each cell is an integer that says which cluster it was
  # just sampled from.  I think that we need to turn the NAs into 0s to make that happen...

  # here is a quick function I need
  sample_from_brick_probs <- function(z) {
    if (any(is.na(z))) return(NA)
    sample.int(n = length(z), size = 1, replace = FALSE, prob = z)
  }

  # then we get the group (color, ultimately) for each cell in a single raster
  pixel_ints <- calc(TRB_rescaled, sample_from_brick_probs)

  # Now we make the rgb raster...

  # here is a matrix in which the columns are rgb values for the colors we want.
  rgb_mat <- col2rgb(cols[names(TRB)])

  # now, we want to add an alpha channel on there too...
  rgba_mat <- rbind(rgb_mat, alpha = floor(alpha * 255))

  # now, a function to pick out columns of that.
  pick_col <- function(x) {
    rgba_mat[, as.integer(x)]
  }

  # this is a little silly that you can't have brick of height one and return a length>1 vector...
  # so we make a brick of two layers and then keep just the first four layers of the result
  silly <- brick(pixel_ints, pixel_ints)
  overdone <- calc(silly, fun = pick_col)
  rgba_rast <- raster::subset(overdone, 1:4)
  names(rgba_rast) <- rownames(rgba_mat)

  if(!is.null(alpha_scale)) {
    mins <- cellStats(scal, min, na.rm = TRUE)
    a <- scale_min_to
    s <- sum(scal_z)
    m <- sum(keep_cells * mins) # raster of the min values in each focal group

    bot_scaled <- ((a + m * (s - m)) / (m + a * (1 - s))) ^ alpha_exp


    rgba_rast[[4]] <-  floor(bot_scaled * alpha_scale * 255)
    rgba_rast[[4]][rgba_rast[[4]] > alpha_chop_max] <- alpha_chop_max
    names(rgba_rast)[4] <- "alpha"
  }

  rgba_rast
}
