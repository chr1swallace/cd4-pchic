library(randomFunctions)
library(devtools)
library(data.table)
library(GenomicRanges)
library(magrittr)

args <- list(d=file.path("/scratch/cew54/FM-GUESS","10p-30686933-30807686"))
source("/home/ob219/git/cd4chic/GUESSFM/R/GUESS-common.R")
CD4CHIC.OUT<-'/home/ob219/scratch/DATA/JAVIERRE_ICHIP/'
## aggregate regions
fs <- list.files("/home/ob219/scratch/DATA/JAVIERRE_ICHIP/RDATA/",full=TRUE,pattern='javierre_tnact')
fs<-fs[grep("gfm",fs,invert=TRUE)]
for(f in fs)
    print(load(f))
## cs.gr = coding snps
## frags.gr = HindIII frags linked to genes as either promoter or interacting

## for this version, we will exclude promoters, see v2 for including them.

## get the GRangesList object that has the intervals in

load("/home/ob219/scratch/DATA/JAVIERRE_ICHIP/RDATA/GUESSFM_intervals_tnact.RData")

logsum <- function(x) {
    my.max <- max(x) ##take out the maximum value in log form)
    my.res <- my.max + log(sum(exp(x - my.max )))
    return(my.res)
}

MPPI.overlap <- function(SM2, in.gr, out.gr) {
    ## extract PP from SM2 over SNPs in in.gr regions, excluding any in out.gr regions
    if(!length(in.gr))
        return(NULL)
    use.gr <- setdiff(in.gr,out.gr)
    if(!length(use.gr))
        return(NULL)
    o <- findOverlaps(use.gr,snps)
    isnps <- subjectHits(o)
    snps.overlap <- snps[isnps,]$snp.name
    lapply(SM2, function(sm) {
        wh <- which(sapply(sm@model.snps, function(ss) any(ss %in% snps.overlap)))
        if(!length(wh))
            return(0)
        logsum(sm@models[wh,"lPP"]) %>% exp()
    }) %>% unlist()
}

splitgr <- function(gr,snps) {
    o <- findOverlaps(gr,snps)
    uq <- unique(queryHits(o)) 
    tmp <- gr[uq]
    with(tmp, split(uq, list(ensg, ld.id)))
}
limitgr <- function(gr,snps,keep=c("ensg","ld.id")) {
    o <- findOverlaps(gr,snps)
    uq <- unique(queryHits(o))
    mcols(gr) <- mcols(gr)[,keep]
    gr[uq]
}

myc <- function(...) {
    l <- list(...)
    l <- l[ sapply(l,length)>0 ]
    do.call("c",l)
}

## TODO: get interaction scores only

## RUN THINGS
library(parallel)


#agg.gr <-GUESS.ichip.grl[['promoter_only']]
for(cat in names(GUESS.ichip.grl)){
	if(cat!='overall'){
		next
	}
	message(paste("Processing",cat))
	agg.gr<-GUESS.ichip.grl[[cat]]
	#agg.gr <- frags.gr[ frags.gr$type=="interaction", ]
	
	## setkey(agg,ensg,ld.id)
	## agg.gr <- with(agg, GRanges(Rle(seqnames),
	##                             ranges=IRanges(start=start,width=width)))
	## mcols(agg.gr) <- agg[,.(ensg,bait,frag.id,ld.id)]
	
	## we will NOT aggregate over frag.id, within .(ensg,ld.id)
	## agg2 <- copy(agg)
	## setkey(agg2,ensg,ld.id)
	## agg2 <- unique(agg2[,.(ensg,ld.id)])
	agg2 <- as.data.table(mcols(agg.gr))
	agg2 <- unique(agg2,by=c("ensg","ld.id"))
	
	DISEASES <- c("ICOELIAC","COELIAC","RA","T1D","GRAVES","JIA","MS", "RA.UK","RA.US","PBC")
	tmp <- matrix(0,nrow(agg2),length(DISEASES),dimnames=list(NULL,DISEASES))
	agg2 <- cbind(agg2,as.data.table(tmp))
	#mcols(agg.gr) <- cbind(mcols(agg.gr), as.data.frame(tmp))
	
	
	dirs <- list.files("/scratch/cew54/FM-GUESS",full=TRUE)
	dirs <- dirs[dir.exists(dirs)]
	sd<-sapply(dirs,function(x){
		if(file.exists(file.path(x,'skip')))
			return(FALSE)
		return(TRUE)
	})
	dirs<-dirs[sd]
	
	
	
	
	for(dd in dirs) {    
	    chr <- sub("[pq].*","",basename(dd))
	    f.geno <- paste0("/scratch/cew54/FM-impute_qc/imputed-",basename(dd),".RData")
    	    f.sm <- file.path(dd, "snpmod-fixed2.RData")
	    if(!file.exists(f.sm))
	        next
	    message("\n",dd)
	    
	## summary of GUESSFM
	    (load(f.sm))
	    (load(f.geno))
	    snps <- with(DATA@snps, GRanges(Rle(rep(chr,nrow(DATA@snps))),
	                                    ranges=IRanges(start=position,width=pmax(nchar(A1),nchar(A2)))))
	    mcols(snps) <- data.frame(snp.name=rownames(DATA@snps),stringsAsFactors=FALSE)
	
	    ## genes/ld.ids overlapping
	    regions <- limitgr(agg.gr,snps)
	    csnps <- limitgr(cs.gr, snps)
	    
	    ## regions.ensg <- sub("\\..*","",names(regions)) #with(tmp, split(as.character(ensg), list(ensg, ld.id))) %>% lapply(., unique)
	    ## regions.ld <- sub(".*\\.","",names(regions)) #with(tmp, split(ld.id, list(ensg, ld.id))) %>% lapply(., unique)
	    ## csnps.ensg <- sub("\\..*","",names(csnps))
	    ## csnps.ld <- sub(".*\\.","",names(csnps))
	
	    all.ensg <- c(regions$ensg, csnps$ensg) %>% unique()
	    all.ld <- c(regions$ld.id, csnps$ld.id) %>% unique()
	    all.ensg.ld <- c(paste(regions$ensg, regions$ld.id, sep="%"),
	                     paste(csnps$ensg, csnps$ld.id, sep="%")) %>% unique()
	    
	    ## check we haven't done this
	    existing <- agg2[ensg %in% all.ensg & ld.id %in% all.ld,
	                     names(SM2),
	                     with=FALSE] %>% as.matrix()
	    if(any(existing>0))
	        warning("revisiting some gene/ld block pair",dd)
	
	    ## aggregate over models
	    message(length(all.ensg.ld))
	    results <- mclapply(all.ensg.ld, function(i) {
	        g <- sub("%.*","",i)
	        l <- sub(".*%","",i)
	        message(l,"\t",g)
	        MPPI.overlap(SM2,
	                     in.gr = myc(regions[ regions$ensg==g & regions$ld.id==l ],
	                                 csnps[ csnps$ensg==g & csnps$ld.id==l ]),
	                     out.gr = csnps[csnps$ensg!=g])
	    }, mc.cores=24)
	
	    ## store results
	    for(i in seq_along(all.ensg.ld)) {
		print(i)
	        g <- sub("%.*","",all.ensg.ld[i])
	        l <- sub(".*%","",all.ensg.ld[i])
		if(!is.null(names(results[[i]]))){
	        	agg2[ensg==g & ld.id==l,
	             names(results[[i]]):=as.list(results[[i]])]
		}
	    }
	
	    if(any(existing>0)) {
	        revised <- agg2[ensg %in% all.ensg & ld.id %in% all.ld,
	                        names(SM2),
	                        with=FALSE] %>% as.matrix()
	        revised <- pmax(revised,existing) %>% as.data.frame() %>% as.list()
	        agg2[ensg %in% all.ensg & ld.id %in% all.ld,
	             names(revised):=revised]    
	    }
	}
	fname<-paste0('MPPI-',cat,'.RData')
	save(agg2,file=file.path(CD4CHIC.OUT, "out/GUESSFM_JUNE_13",fname))
}

