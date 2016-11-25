genoscapeRtools
================
25 November, 2016

-   [Introduction](#introduction)
-   [Data conversion](#data-conversion)
-   [Running admixture with a single function.](#running-admixture-with-a-single-function.)

<!-- README.md is generated from README.Rmd. Please edit that file -->
Introduction
------------

The `genoscapeRtools` has a few functions that simplify the process of running the program `admixture` multiple times and multiple K values. This is designed to run on a Mac with a decent Unix-like system. Hasn't been tested elsewhere.

Data conversion
---------------

The first hurdle is to get your data into a proper PLINK "bed" file that `admixture` can read. This is pretty easy from the VCF file, but we are going to do it straight from the object that `read_012` creates.

So, read that thing in:

``` r
library(genoscapeRtools)

wifl_clean <- read_012("../cleaned_indv175_pos105000")
```

Then, if you want to make a plink bed file out of that you can do it like this:

``` r
convert_012_to_bed(wifl_clean, chromo_override = TRUE, prefix = "plinky")
#> converting 012 genotypes to plink ped format
#> writing plink ped file to plinky.ped
#> writing plink map file to plinky.map
#> using plink to create files plinky.{bed,bim,fam}
```

That creates all these files:

``` r
dir(pattern = "plinky.*")
#> [1] "plinky.bed" "plinky.bim" "plinky.fam" "plinky.log" "plinky.map"
#> [6] "plinky.ped"
```

And now you can pass the `plinky.bed` file to our function that runs admixture a lot of times. That function will look for the .fam and the .bim files in the same place, so they have to be there too.

Running admixture with a single function.
-----------------------------------------

This function sets up a directory and runs admixture in parallel in all of them. We use the defaults, which puts all the results in a directory called "admixture\_runs" in the current working directory.

We are going to just do one rep and only at K=3 and K=4

``` r
run_admixture(bed = "plinky.bed", Reps = 1, Kvals = 3:4)
```

The above can take a long time and it is hard to kill the process once you start it. So, just let 'er go.

To see how things are coming along you can look at the files inside the `r_x_K_y` directories in the `admixture_runs` directories where all the action is happening.

When it is all done, the output is in the directories like this:

    2016-11-25 06:34 /tutorials/--% (master) ls admixture_runs/*
    admixture_runs/data:
    README.txt  input.bed   input.bim   input.fam

    admixture_runs/r_1_K_3:
    admixture.stderr  admixture.stdout  input.3.P         input.3.Q

    admixture_runs/r_1_K_4:
    admixture.stderr  admixture.stdout  input.4.P         input.4.Q

So, all that remains is a simple function to read all that output into a nice tidy data frame.

And, I need to add a cross-validation flag if we want to do cross-validation.
