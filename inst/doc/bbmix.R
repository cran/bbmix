## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----call_gt------------------------------------------------------------------
## Retrive input files for running call_gt
counts_f <- system.file("extdata/input", "NA12878.chr22.Q20.allelicCounts.txt",
		package = "bbmix",
		mustWork = TRUE)

legend <- system.file("extdata/input", "1000GP_Phase3_chr22.legend",
       package = "bbmix", mustWork = TRUE)

fisher_f <- system.file("extdata/input", "chr22.FS.Q20.alleleCounts.txt",
	 package = "bbmix", mustWork = TRUE)

cluster_f <- system.file("extdata/input", "fSNPs_22_RP_maf0_01_cluster3window35.txt",
	  package = "bbmix", mustWork = TRUE)

library(bbmix)

## Choose your output file:

out <- paste0(tempdir() , "/NA12878.chrom22.gt.txt")

## Run call_gt:
call_gt(allele_counts_f = counts_f,
	legend_f = legend,
	fisher_f = fisher_f,
	cluster_f = cluster_f,
	out = out)
unlink(out)

## ----get_files----------------------------------------------------------------
## Retrive input files for running call_gt
counts_f <- system.file("extdata/input", "NA12878.chr22.Q20.allelicCounts.txt",
		package = "bbmix",
		mustWork = TRUE)


## ----fit_bb, eval=F-----------------------------------------------------------
#  ## Choose your output directory, in this examle we use a temporary directory, also we use N=10 for the model to run quickly but this should be 1000 for model convergence (function default)
#  
#  out <- "path/to/output_dir"
#  out <- tempdir()
#  
#  ## Run fit_bb:
#  fit_bb(counts_f = counts_f,
#  	N=10,
#  	out = out, mc.cores=1)
#  unlink(out)

## ---- eval=F------------------------------------------------------------------
#  ## Run call_gt:
#  call_gt(allele_counts_f = counts_f,
#  	stan_f = output_from_fit_bb,
#  	legend_f = legend,
#  	fisher_f = fisher_f,
#  	cluster_f = cluster_f,
#  	out = out)

## ----poolreads----------------------------------------------------------------
## Run poolreads:

counts_f <- system.file("extdata/input", "NA12878.chr22.Q20.allelicCounts.txt",
		package = "bbmix",
		mustWork = TRUE)

## In this example we only use one file
out <- tempfile()

poolreads(files=counts_f,
	N=10,
	d=10,
	out = out)
	
unlink(out)

## ----inspect_genotype_calls---------------------------------------------------
## Inspecting output files
gt <- data.table::fread(system.file("extdata/output", "gt.NA12878.chr22.txt",
      						      package = "bbmix",
						      mustWork = TRUE))

head(gt)


## ----qc_genotype--------------------------------------------------------------

library(bbmix)
## Extracting genotype file 
gt_f <- system.file("extdata/output", "gt.NA12878.chr22.txt",
      						      package = "bbmix",
						      mustWork = TRUE)	
out <- tempfile()

## Running function
ex_alt_hom(gt_f, out)

unlink(out)

