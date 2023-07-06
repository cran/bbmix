rstan::rstan_options(auto_write = TRUE)




##' Fit beta binomial distribution to allelic counts for homozygous reference, heterozygous, homozygous alternative
##' 
##' @param counts_f file name with allele counts for SNPs
##' @param depth  depth cut-off to use to select SNPs to fit  distributions
##' @param N number of SNPs to use for fitting
##' @param k number of components for mixture model, defaults to 3
##' @param alpha_p alpha parameter for the k components of alpha parameter 
##' @param beta_p beta paramenter for the k components of Beta parameter
##' @param prefix charcter with prefix to add for saving files, defaults to NULL
##' @param mc.cores number of cores to use, defaults to parallel detected cores
##' @param out character with dir name to save output
##' @import data.table Rcpp
##' @importFrom rstan sampling
##' @importFrom utils write.table
##' @export
##' @return saves stan object to file
##' @examples
##'
##' \dontrun{
##' ## Retrive input files for running call_gt
##' counts_f <- system.file("extdata/input", "NA12878.chr22.Q20.allelicCounts.txt",
##' package = "bbmix",
##' mustWork = TRUE)
##'
##' out <- tempdir()
##' fit_bb(counts_f = counts_f, N=10,
##' out = out, mc.cores=1)
##' unlink(out)
##' }
##' 
fit_bb <- function(counts_f, depth=10, N=1000, prefix=NULL, k=3, alpha_p=c(1,10,499), beta_p=c(499, 10, 1), out, mc.cores=NULL) {

    old <- options()
    on.exit(options(old))
    
    if(is.null(mc.cores)){
        options(mc.cores = parallel::detectCores())
    }
    

    Total <- NREF <- NALT <- mu  <- NULL
    
    if(!file.exists(counts_f)) stop(paste("invalid file name for  ", counts_f))
    if(!dir.exists(out)) stop(paste("invalid directory name for  ", out))
    
    
    theta_p <- rep(1,k)
    counts <- fread(counts_f)
    
    ## Apply depth filter to counts and select N rows at random
    counts <- counts[, Total := NREF+NALT][Total >= depth]

    if(N > nrow(counts)) stop(paste("The number of SNPs for model fitting", N, "has to be smaller that the total number of reads in the file: ", nrow(counts)))
    
    counts <- counts[sample(1:nrow(counts), N),]
    
    
    mod_stan <- stanmodels$BetaBinomMix
    
    inp <- lapply(k, function(i) list(N=nrow(counts), K= i, n=counts[['NALT']], m=counts[['Total']], alpha_p=alpha_p, beta_p=beta_p, theta_p=theta_p))

    message("Running stan model")

    ## set initial values for ordered mu  to be between 0 and 1 and linked lambda
    init_mu_lambda <- lapply(1:4, function(i) {
        tmp <- data.table(mu=rbeta(k, alpha_p, beta_p),lambda=rgamma(k, alpha_p + beta_p, 1))
        setkey(tmp,mu)       
        return(list(mu=tmp$mu, lambda=tmp$lambda))
    })


    stan_fit <- lapply(inp, function(i) rstan::sampling(mod_stan,data=i, cores=mc.cores, refresh=0 , init=init_mu_lambda))
    

    if(is.null(prefix)) {
        saveRDS(stan_fit, paste0(out, "/stan_beta_binom_fit.rds"))
        

    } else {
        saveRDS(stan_fit, paste0(out, "/", prefix, "_stan_beta_binom_fit.rds"))
    
    }
    
}

##' Call genotypes using beta binomial after model training
##' 
##' @param allele_counts_f vector with file names with allele counts for SNPs
##' @param depth min read count to call variant
##' @param stan_f full name to stan object with model fit to extract mean of parameters. Defaults to the model trained with genome in a bottle reads. Otherwise this object can be generated with fit_bb function.
##' @param legend_f full name for file with SNP info to get allele frequency for prior
##' @param pop population to select AF for GT prior, defaults to EUR
##' @param prob cut-off for making hard calls, defaults to 0.99
##' @param fisher_f file with Fisher test to detect strand bias
##' @param fisher cut_off for Fisher test to detect strand bias
##' @param cluster_f file with info about SNP clusters
##' @param out character with file name to save genotype output
##' @import data.table R.utils
##' @importFrom utils write.table
##' @export
##' @return data table with genotype probabilities
##' @examples
##' ## Retrive input files for running call_gt
##' counts_f <- system.file("extdata/input", "NA12878.chr22.Q20.allelicCounts.txt",
##' package = "bbmix",
##' mustWork = TRUE)
##' 
##' legend <- system.file("extdata/input", "1000GP_Phase3_chr22.legend",
##' package = "bbmix", mustWork = TRUE)
##' 
##' fisher_f <- system.file("extdata/input", "chr22.FS.Q20.alleleCounts.txt",
##' package = "bbmix", mustWork = TRUE)
##' 
##' cluster_f <- system.file("extdata/input", "fSNPs_22_RP_maf0_01_cluster3window35.txt",
##' package = "bbmix", mustWork = TRUE)
##'
##' out <- paste0(tempdir() , "/NA12878.chrom22.gt.txt")
##'
##' ## Run call_gt:
##' call_gt(allele_counts_f = counts_f,
##' legend_f = legend,
##' fisher_f = fisher_f,
##' cluster_f = cluster_f,
##' out = out)
##'
##' unlink(out)
##' 
call_gt <- function(allele_counts_f, depth=10, stan_f=NULL, legend_f, pop='EUR', prob=0.99, fisher_f=NULL, fisher=30, cluster_f = NULL, out){

    total <- NREF <- NALT <- p0 <- p1 <- p2 <- GT <- CHROM <- SNPcluster <- FS <- NULL
    ## check inputs:
    files <- c(allele_counts_f, legend_f)
    if(!is.null(stan_f)){
        files <- c(files, stan_f)
    }

    if(!is.null(cluster_f)) files <- c(files, cluster_f)

    w <- !file.exists(files)
    if(any(w)) stop(paste("invalid file names ", paste(files[w], collapse= ", ")))

    ## Get model parameters
    
    
    if(is.null(stan_f)) {              
        stan_samples <- gbfit             #default internal object
    } else {
        pars <- c("mu", "lambda")
        stan_fit <- readRDS(stan_f)
        stan_samples <- rstan::extract(stan_fit[[1]], pars)
    }
    
   
    ## Get counts
    counts <- fread(allele_counts_f)
    counts[ , total := NREF + NALT]
    counts <- counts[total >= depth,]

    ## Get allele frequency to calculate P(G=g | D)
    
    snps <- fread(legend_f)
    

    ## Make sure to exclude SNPs based on eaf
    counts <- merge(counts, snps[, c("position", paste0('a', 0:1), pop), with=F], by.x=c("POS", "REF", "ALT"), by.y=c("position",  paste0('a', 0:1)))

    setcolorder(counts, c("CHROM", "POS", "REF", "ALT"))

    ## Get mean posterior for each genotype, mean expected value and sd expected value for each fSNP
    post <- rbindlist(lapply(1:nrow(counts), function(i) gt_help(stan_samples, pop, counts[i,])))

    ## write.table(cbind(counts, post), file=out, row.names=F)

    gt <- cbind(counts, post)

    ## Make hard calls based on prob

    gt[p0 >prob | p1>prob | p2 >prob,  GT := names(.SD)[max.col(.SD)], .SDcols = paste0("p", 0:2)][, GT := as.integer(gsub("p", "", GT))]

    ## perform QC

    if(!is.null(cluster_f)){            #remove clustered SNPs

        clusters <- fread(cluster_f)

        clusters[, CHROM := as.integer(gsub(".*fSNPs_([0-9]+)_RP.*", "\\1", cluster_f))]

        gt <- merge(gt, clusters, by=c("CHROM", "POS", "REF", "ALT"))

        gt <- gt[SNPcluster == FALSE,]

    }

    if(!is.null(fisher_f)){

        fish <- fread(fisher_f)

        gt <- merge(gt, fish[, c("CHROM", "POS", "REF", "ALT", "FS")], by=c("CHROM", "POS", "REF", "ALT"))
        gt <- gt[FS  < fisher,]

    }

    write.table(gt, file=out, row.names=F)

}

##' call gt helper, get posterior mean, expected gt and sd expected gt across all samples
##'
##' @param stan_samples matrix with samples extracted from stan fit object, params mu and lambda
##' @param pop population to select AF for GT prior, defaults to EUR
##' @param data data table 1 row with counts and EAF to apply model
##' @return
##' gt_help()
gt_help <- function(stan_samples, pop, data){
    
    ## get likelihood in matrix form
    d <- Reduce(cbind, lapply(1:3, function(i) rmutil::dbetabinom(y=data[['NALT']], size=data[['total']], m=stan_samples$mu[,i], s=stan_samples$lambda[,i])))

    ## Apply to each column prior (p(G==g)) based on EAF
    post <- d %*% diag(c((1 - data[[pop]])^2, 2*data[[pop]]*(1-data[[pop]]), data[[pop]]^2) )

    ## Normalise each posterior by the sum

    sumP <- rowSums(post)
    postN <- post * 1/sumP

    ## Calculate mean for each posterior mean

    postN_mean <- colMeans(postN)

    ## calculate expected value

    expV <- rowSums(postN %*% diag(c(0,1,2)))

    ## prepare data table to return
    dt <- as.data.table(matrix(c(postN_mean, mean(expV), sd(expV)), nrow=1, dimnames=list(NULL, c(paste0("p", 0:2), "expected_GT", "expected_sd"))))

    return(dt)

}


##' call gt helper, calculate mean dbetabinom from all posterior samples
##'
##' @param n counts alt allele
##' @param m total counts
##' @param mu vector with posterior draws for mu param
##' @param lambda vector with posterior draws for lambda param
##' @return mean of dbetabinom
##' 
call_help <- function(n, m, mu, lambda){

    return(mean(rmutil::dbetabinom(n, m, mu, lambda)))
    
}

##' Pool randomly selected reads from different files
##'
##' @param files names for files to extract reads
##' @param N number of reads to extract
##' @param d depth for reads
##' @param out file name to save reads
##' @import data.table
##' @importFrom utils write.table
##' @export
##' @return save files
##' @examples
##' counts_f <- system.file("extdata/input", "NA12878.chr22.Q20.allelicCounts.txt",
##' package = "bbmix",
##' mustWork = TRUE)
##' 
##' ## In this example we only use one file and we take a pool of 10 reads
##' 
##' out <- tempfile()
##' 
##' poolreads(files=counts_f,
##' N=10,
##' d=10,
##' out = out)
##'
##' unlink(out)
##' 
poolreads <- function(files, N=1000, d=10,  out){

    NREF <- NALT <- NULL
    
    n=1
    for( f in files){
        dt <- fread(f)
        dt <- dt[NREF + NALT >= d,]
        dt <- dt[sample(1:nrow(dt), N),]
        if(n == 1) {
            write.table(dt, out, row.names=F)

        } else {

             write.table(dt, out, row.names=F, append=T, col.names=F)
        }
        
        
        n=2

    }

}

##' Exclude fSNPs with no alternative allele in any sample. Also exclude fSNPs if all samples are hom.
##'
##' @param gt_f character vector with file names with genotype calls per sample
##' @param out file name to save output
##' @import data.table
##' @importFrom utils write.table
##' @export
##' @return save file
##' @examples
##' 
##' gt_f <- system.file("extdata/output", "gt.NA12878.chr22.txt",
##' package = "bbmix",
##' mustWork = TRUE)
##' out <- tempfile()
##'
##' ## Running function
##' ex_alt_hom(gt_f, out)
##'
##' unlink(out)
##' 
ex_alt_hom <- function(gt_f, out){

    id <- CHROM <- POS <- REF <- ALT <- p0  <- p1  <- p2  <- hard <- NULL

    gt <- lapply(gt_f, fread)

    names(gt) <- lapply(gt_f, function(i) gsub('([A-Z]+[0-9]+).chrom.*', '\\1', basename(i)))

    gt_long <- rbindlist(gt, idcol='sample')

    gt_long[, id:= paste(CHROM, POS, REF, ALT, sep="_")]
    
    gt_hard <- as.matrix(gt_long[ , list(p0, p1, p2)])

    ## Get hard calls, element with max prob minus 1 to make scale 0,1,2
    max <- apply(gt_hard, 1, function(i) which(i==max(i)) -1 )

    gt_long[, hard := max]

    ## Get ids for homs across all samples

    wide <- dcast(gt_long, id ~ sample, value.var='hard') 

    comp <- wide[complete.cases(wide), ]

    ## find hom in all samples and get id col
    rem <- sort(unlist(lapply(c(0,2), function(i) comp[which(rowSums(comp[, !"id"]) == i*(ncol(comp) - 1)), id] )))
    
    ## exclude from gt_long

    gt_clean <- gt_long[!id %in% rem,]

    ## exclude fSNPs if NA or 0 across all samples

    incomp <- wide[!complete.cases(wide), ]
    rem_id <- incomp[which(rowSums(incomp[,!"id"], na.rm=T) == 0), id]

    gt_clean <- gt_clean[!id %in% rem_id,][, id:=NULL]

    write.table(gt_clean, out, row.names=F)
    
}



        

    
    
    
    
    

    

    
    
       
        

    

    
    
