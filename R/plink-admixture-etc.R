

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

    map <- cbind(chr, locs, 0, bp)
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
#' @param Reps
#' @param Kvals
#' @param path  the path that gets you to the directory you want to put all the results in.
#' @param outdir the single name of the output directory you want
#' @export
run_admixture <- function(bed, Reps = 2, Kvals = 2:5, path = ".", outdir = "admixture_runs") {
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
  })

}
