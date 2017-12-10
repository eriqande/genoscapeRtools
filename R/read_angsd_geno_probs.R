

#' read an angsd binary genotype probability file
#'
#' This reads the probabilities into an array that is 3 x num_ind x num_snps.
#' Things get returned in a list that has includes a vector of sample names
#' and a vector of locus names.
#' @param gz_geno_file  path to the genotype file.  This is the binary file of
#' doubles that holds the posterior probabilities for each individuals genotype
#' at each locus, and it has to be gzipped.  ie., results.genos.gz
#' @param bamlist path to the file that was the bamlist.  It is assumed that
#' the names of the bam files are the names of the samples.  All the path information
#' is removed (if it exists) and .bam is removed from the end (if it is there.)
#' @param mafs_file the mafs file that comes out of ANGSD.  This is used to get
#' information (especially the chrom and position of each variant). It can be gzipped if desired.
#' @details The geno probabilities come out in an array indexed by [genotype, sample, snp].
#' @export
read_angsd_geno_probs <- function(gz_geno_file, bamlist, mafs_file) {

  # read the bamlist and learn how many individuals there are
  samples <- scan(bamlist, what = "character") %>%
    stringr::str_replace("^.*/", "") %>%
    stringr::str_replace("\\.bam$", "")

  n_samples <- length(samples)

  # read the mafs file
  mafs_tibble <- readr::read_tsv(mafs_file)

  # now we can predict how many doubles there should be in the genos file
  n_doubles <- n_samples * nrow(mafs_tibble) * 3

  # now read the genos.  We prepare for 10% more space than needed, and then
  # we check to make sure that we read the right number of things.
  dvec <- readBin(gzcon(gzfile(gz_geno_file)), what = "double", n = floor(n_doubles * 1.1))

  if (length(dvec) != n_doubles) {
    stop("Was expecting ", n_doubles, " entries in the gz_geno_file. Have read in ", length(dvec),
         ". This indicates a discrepancy between the number of samples in bamlist or the number of loci in the mafs_file,",
         " relative to the number of genos with posteriors in the gz_geno_file.")
  }

  darray <- array(dvec, dim = c(3, n_samples, nrow(mafs_tibble)))

  list(samples = samples, snps = mafs_tibble, probs = darray)
}
