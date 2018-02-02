

# functions for interfacing with PLINK and admixture


#' create a PLINK bed file from and 012 matrix
#'
#' This clearly does not preserve Ref and Alt alleles, etc., but can
#' be a quick way to get to a bed file.  It also allows the user to
#' override the chromosome names so that they play nicely with ADMIXTURE.
#' This function expects that plink is on the system PATH somewhere.
#'
#' This function assumes everyone is unrelated and nobody has a phenotype
#' to be coded in plink.
#'
#' It also currently is hardwired to deal with contig names the way they come
#' out of our pipeline and get catenated by \code{\link{read_012}}.  That is,
#' when it sees: \code{scaffold2|size4533513--4321992} it knows that the contig
#' comes before the -- and the position in bp comes after that.
#' @param W the matrix that results from \code{\link{read_012}}.
#' @param chromo_override If TRUE, this sets the chromosome names to 1
#' so that ADMIXTURE will work.
#' @param prefix Prefix for the output files.
#' @export
convert_012_to_bed <- function(W, chromo_override = FALSE, prefix = "plink") {

    # get individual IDs and locus names
    ids <- rownames(W)
    locs <- colnames(W)

    # get the ID columns of the ped file
    pedid <- cbind(ids, ids, 0, 0, 0, 0)


    # make a bunch of genotypes 0 => "1 1"; 1 => "1 2"; 2 => "2 2", and -1 => "0 0"
    message("converting 012 genotypes to plink ped format")
    g <- W
    mode(g) <- "character"
    g[g == "-1"] = "0\t0"
    g[g == "0"] = "1\t1"
    g[g == "1"] = "1\t2"
    g[g == "2"] = "2\t2"
    dimnames(g) <- NULL

    # write the file
    message(paste0("writing plink ped file to ", prefix, ".ped"))
    write.table(cbind(pedid, g), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, file = paste0(prefix, ".ped"))

    # now put together the map file
    tmp <- stringr::str_split_fixed(locs, "--", 2)
    chr <- tmp[,1]
    bp <- tmp[,2]
    if(chromo_override == TRUE) {
      chr <- 1
    }

    map <- cbind(chr, locs, 0, 1:length(locs))  # put the position as just the index to keep everything in original order
    colnames(map) <- NULL
    message(paste0("writing plink map file to ", prefix, ".map"))
    write.table(map, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, file = paste0(prefix, ".map"))

    # now use plink to convert that to a bedfile
    message(paste0("using plink to create files ", prefix, ".{bed,bim,fam}"))
    system(paste0("plink --file ", prefix, " --aec --out ", prefix))
}




#' run admixture Reps times at each of K values
#'
#' Just a wrapper to mclapply over all this.
#' @param bed  the bed file to use. (With the .bed extension.  Example = "plink.bed".)  The .bim file has to be with it
#' in the same directory as well as the .fam file.
#'  These three files will both get
#' copied to "input.bed" and "input.bim" in the "data" directory in the outdir
#' @param Reps number or reps to do
#' @param Kvals vector of K-values to use
#' @param path  the path that gets you to the directory you want to put all the results in.
#' @param outdir the single name of the output directory you want
#' @param num_cores  the number of cores to parallelize the admixture runs over.  By default it will
#' try to detect the number of cores and then use 1 less than that.  (or 1, whichever is larger). But
#' you can set it to whatever...
#' @export
run_admixture <- function(bed, Reps = 2, Kvals = 2:5, path = ".", outdir = "admixture_runs", num_cores = max(1, parallel::detectCores() - 1)) {
  OUT <- file.path(path, outdir)
  DAT <- file.path(OUT, "data")
  if(dir.exists(OUT)) stop(paste0("Directory ", OUT, " already exists! Remove it if you want to proceed"))

  dir.create(OUT, recursive = TRUE)  # create the output directory
  dir.create(DAT)
  file.copy(bed, file.path(DAT, "input.bed"))
  bim <- stringr::str_replace(bed, "\\.bed$", ".bim")
  file.copy(bim, file.path(DAT, "input.bim"))
  fam <- stringr::str_replace(bed, "\\.bed$", ".fam")
  file.copy(fam, file.path(DAT, "input.fam"))
  cat(paste0("input files copied here from ", bed, "\n"), file = file.path(DAT, "README.txt"))

  # now make a matrix of calculation directories, rep numbers, and K values
  Rs <- rep(1:Reps, each = length(Kvals))
  Ks <- rep(Kvals, Reps)
  dirs <- cbind(file.path(OUT, paste0("r_", Rs, "_K_", Ks)), Rs, Ks)

  # create the directories to the calculations
  dump  <- lapply(dirs[, 1], function(x) dir.create(x))

  # run admixture with a different seed each time
  parallel::mclapply(1:nrow(dirs), function(i) {
    system(paste0("cd ", dirs[i, 1], "; admixture  -s ", i * 1234, " ../data/input.bed ", dirs[i, 3], " > admixture.stdout 2> admixture.stderr "))
  }, mc.cores = num_cores
  )

}




#' read in all the AMIXTURE files that are in the r_X_K_Y directories produced by run_admixture()
#'
#' Reads it into a tidy data frame
#' @param path  the directory that holds all the r_X_K_y directories.
#' @param which_dat  By default, this gets the Q values.  If you want to get
#' the allele frequencies, set this to "P"
#' @export
slurp_admixture <- function(path = "admixture_runs", which_dat = "Q") {
  stopifnot(which_dat %in% c("P", "Q"))
  dirs <- dir(path = path)
  dirs <- dirs[stringr::str_detect(dirs, "^r_[0-9]+_K_[0-9]+")]  # get the results directories

  # read the loci in
  bimfile <- file.path(path, "data", "input.bim")
  bim <- readr::read_delim(bimfile, delim = "\t", col_names = FALSE)

  locs <- bim$X2  # vector of locus names

  famfile <- file.path(path, "data", "input.fam")
  fam <- readr::read_delim(famfile, delim = "\t", col_names = FALSE)

  ids <- fam$X1

  Ks <- stringr::str_split_fixed(dirs, "_", n = 4)[,4]
  Rs <- stringr::str_split_fixed(dirs, "_", n = 4)[,2]

  ret <- lapply(1:length(dirs), function(i) {
    dir <- dirs[i]
    K <- Ks[i]
    rep <- Rs[i]

    if(which_dat == "Q") {
      dat <- readr::read_delim(file.path(path, dir, paste0("input.", K, ".Q")), delim = " ", col_names = FALSE) %>%
        setNames(paste0("c_", 1:K)) %>%
        mutate(ID = ids, K = K, rep = rep) %>%
        select(ID, K, rep, everything()) %>%
        tidyr::gather(key = cluster, value = "Q", -ID, -K, -rep)
    } else {
      dat <- readr::read_delim(file.path(path, dir, paste0("input.", K, ".P")), delim = " ", col_names = FALSE) %>%
        setNames(paste0("c_", 1:K)) %>%
        mutate(pos = locs, K = K, rep = rep) %>%
        select(pos, K, rep, everything()) %>%
        tidyr::gather(key = cluster, value = "P", -pos, -K, -rep)
    }
    dat
  }) %>%
    dplyr::bind_rows()

}



#' ggplot the Qs from multiple runs of ADMIXTURE
#'
#' This just plots everything onto one page and returns it as a ggplot object.
#' @param Qs the data frame that comes from running \code{\link{slurp_admixture}}
#' on the Qs.
#' @param sort_df a data frame that has, at a minimum, a column called "ID" that has the
#' birds you want to include, in the order you want to include them.
#' @export
ggplot_the_Qs <- function(Qs, sort_df = NULL) {

  if(!is.null(sort_df)) {
    bird_levels <- sort_df$ID
    missers <- setdiff(unique(Qs$ID), bird_levels)
    if(length(missers > 0)) {
      warning("Sorting birds with sort_df will remove the following non-ordered individuals: ", paste(missers, collapse = ", "))
    }
    Qs$ID <- factor(Qs$ID, levels = bird_levels)
    Qs <- Qs %>% filter(!is.na(ID))
  }
  # get the colors that distruct uses by default (I got these by pulling them out of the PS code.  You gotta set anything larger than 1 to 1.)
  distruct_colors <- c("#FF994D", "#0099E6", "#E6FF00", "#FF99E6", "#339933", "#800080", "#FF004D", "#00FF00", "#0000FF", "#FF00FF", "#FFE699", "#B34D00", "#00FFFF", "#808000", "#FF9999", "#008080", "#99BF26", "#7326E6", "#26BF99")
  maxK <- max(Qs$K)
  g <- ggplot(Qs, aes(x = ID, y = Q, fill = cluster)) +
    geom_bar(stat = "identity", position = "stack", width = 1.0) +
    scale_fill_manual(values = distruct_colors[1:maxK]) +
    facet_grid(K + rep ~ ., drop=TRUE, space="free", scales="free") +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 0.1))

  g
}


