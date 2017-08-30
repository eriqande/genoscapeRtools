

#' do a PCA using SNPRelate and prepare output for rapid plotting
#'
#' This just does all the legwork to do a PCA using SNPRelate, and then
#' it does what it takes to plot every PC against every other using
#' facet_grid.  It doesn't make the plots, but it prepares the
#' data.
#'
#' This function lets you say which individuals you want to drop which
#' makes it easy to iteratively remove outliers.
#' @param dat012 an 012 matrix like that you get from read_012.
#' missing data should be denoted -1.
#' @param samples_to_remove the names of the samples that should be removed
#' before doing the PCA. This does not check to make sure that the samples
#' requested to be dropped actually exist in dat012.
#' @param num_pcs number of principal components to retain in output.  Default = 6
#' @param sample_groups If not NULL, this should be a tibble with a column "sample"
#' and another column "group", by which the points will be colored.
#' @param plot_facet_free if TRUE then the pca-pairs plot scales are free to encompass less area.
#' @export
#' @examples
#' Just while working on it, for testing:
#'
#' dat012 <- read_012("../amke-popgen/data/rachael-amke-clean_indv175_pos110000", gz = TRUE)
#' sample_groups <- read_csv("../amke-popgen/data/amke_locations_from_mikki.csv")
#' samples_to_remove <- c("1833-14592", "2003-40521", "OCBPC-3", "OCBPC-4")
genoscape_pca <- function(dat012,
                          samples_to_remove = character(0),
                          num_pcs = 6,
                          sample_groups = NULL,
                          plot_facet_free = FALSE) {

  # a little error checking:
  if (!is.null(sample_groups) && !("group" %in% names(sample_groups))) {
    stop("If you supply a sample_groups data frame it must have a column named \"groups\"")
  }

  # first deal with removing individuals and telling the user about it
  indivs <- rownames(dat012)
  droppers <- indivs[indivs %in% samples_to_remove]
  if (length(droppers) > 0) {
    message("Dropping: ", paste(droppers, collapse = ", "))
  }
  keepers <- indivs[!(indivs %in% samples_to_remove)]

  new012 <- dat012[keepers, ]

  # next do the PCA.  Write gds files to a temporary location
  gds_file <- tempfile()

  SNPRelate::snpgdsCreateGeno(gds.fn = gds_file,
                   genmat = new012,
                   sample.id = rownames(new012),
                   snp.id = colnames(new012),
                   snpfirstdim = FALSE)

  gds_conn <- SNPRelate::snpgdsOpen(gds_file)
  sr_pca <- SNPRelate::snpgdsPCA(gds_conn, autosome.only = FALSE)
  SNPRelate::snpgdsClose(gds_conn)

  # now make a short-format tibble of the pca results
  pca_tib <- dplyr::bind_cols(
    tibble::tibble(sample = sr_pca$sample.id),
    tibble::as_tibble(sr_pca$eigenvect[, 1:num_pcs])
  ) %>%
    stats::setNames(c("sample", sprintf("PC-%02d", 1:num_pcs)))

  # and make a long version that has a PCx and a PCy.
  # first make a long tibble
  pca_long <- pca_tib %>%
    tidyr::gather(., key = "PC", "val", -sample)

  # then expand a grid of the possible comparisons (ordered)
  expg <- expand.grid(sample = pca_tib$sample,
                      PCx = sprintf("PC-%02d", 1:6),
                      PCy = sprintf("PC-%02d", 1:6),
                      stringsAsFactors = FALSE) %>%
    tibble::as_tibble()

  # then left join the pca results onto that
  pca_pairs <- dplyr::left_join(expg, pca_long, by = c("sample", "PCx" = "PC")) %>%
    dplyr::rename(val_x = val) %>%
    dplyr::left_join(pca_long, by = c("sample", "PCy" = "PC")) %>%
    dplyr::rename(val_y = val)


  # now, let's associate the meta data to them if it is not null
  if (!is.null(sample_groups)) {
    pca_tib <- dplyr::right_join(sample_groups, pca_tib, by = "sample")
    pca_pairs <- dplyr::left_join(pca_pairs, sample_groups, by = "sample")
  }

  if (!is.null(sample_groups)) {
    g <- ggplot(pca_pairs, aes(x = val_x, y = val_y, fill = group))
  } else {
    g <- ggplot(pca_pairs, aes(x = val_x, y = val_y))
  }

  g2 <- g +
    geom_point(pch = 21, size = 2) +
    scale_fill_discrete(na.value = "white")

  if (plot_facet_free == FALSE) {
    pairs_plot <- g2 +
      facet_grid(PCy ~ PCx)
  } else {
    pairs_plot <- g2 +
      facet_grid(PCy ~ PCx, scales = "free")
  }

  # then return it all
  list(
    pca_df = pca_tib,
    pca_pairs = pca_pairs,
    pairs_plot = pairs_plot,
    dropped_indivs = tibble::tibble(sample = droppers),
    requested_to_drop = tibble::tibble(sample = samples_to_remove)
  )

}
