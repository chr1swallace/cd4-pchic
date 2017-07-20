#!/usr/bin/env Rscript
## follows on from GUESS-qc.R, may need some code from that inserted

library(randomFunctions)
library(annotSnpStats)
library(magrittr)
library(RColorBrewer)
library(reshape)
library(plyr)
library(igraph)
library(ggplot2)
theme_set(theme_bw())
library(gridExtra)
setwd(file.path(CD4CHIC.ROOT,"GUESSFM"))
source("~/DIRS.txt")

## files
args <- getArgs(default=list(d=file.path(ROOT,"FM-GUESS/10p-6030243-6169685"),
                             cond=FALSE))
myfile <- function(stub) {
    if(args$cond)
        stub <- paste0("cond-",stub)
    file.path(plot.dir,stub)
}
## functions and set some parameters
source("R/GUESS-common.R")

## HAPLOTYPES
(load(file.path(args$d,"snpmod.RData"))) # SM2 snpmods
(load(file.path(args$d,"snp-picker-grouped.RData")))
library(snpHaps)
library(snpStatsWriter)
DISEASES <- intersect(c("T1D","ICOELIAC","COELIAC","GRAVES","MS","JIA","RA","RA.UK","RA.US"),
                      names(SM2))

## again <- TRUE; thr <- 0.1
## while(again) {
bm <- best.snps(SM2[DISEASES])
##   thr <- thr/2
##   again <- any(sapply(bm,nrow)==0)
## }

bmsnps <- lapply(bm,rownames) %>% unlist() %>% unique() %>% setdiff(.,"1")
hsnps <- as(gSP2,"groups") %>% tags()
hsnps <- unique(c(bmsnps,hsnps))
lapply(SM2, function(x) x@snps[hsnps,])

if(args$cond) {
    (load(file.path(args$d,"conditional.RData")))
    if(nrow(conditional))
        hsnps <- unique(c(hsnps,as.character(subset(conditional,p<1e-6,select="snp",drop=TRUE))))
}

if(basename(args$d)=="19p-10396336-10628468")
hsnps <- setdiff(hsnps,"rs2278442.10444826.G.A") # only found in RA.US, not believed

## allele frequencies
(load(f.geno))
p <- DATA@samples
if(all(grepl("samples",colnames(p)))) {
  colnames(p) <- sub("samples.","",colnames(p))
  samples(DATA) <- p
}
ph <- unique(DATA@samples[,c("phenotype","country")])
ph <- subset(ph,!is.na(phenotype))
ph <- ph[order(ph$phenotype,ph$country),]
maf <- matrix(0,nrow(ph),length(hsnps), dimnames=list(apply(ph,1,paste,collapse="/"),
                                          hsnps))
for(i in 1:nrow(ph)) {
  maf[i,] <- col.summary(DATA[p$country==ph[i,"country"] & p$phenotype==ph[i,"phenotype"], hsnps])[,"RAF"]
}
write.csv(maf,file=myfile("allele-frequencies.csv"))

## add PCs to samples data.frame
(load(file.path(ROOT,"FM-pca/PCS+.RData")))
rownames(pcs) <- make.names(rownames(pcs))
df.samples <- pcs[rownames(samples(DATA)), ]
df.samples$y <- as.numeric(df.samples$phenotype!="CONTROL")

#' generate haplotypes
X <- as(DATA[,hsnps],"SnpMatrix")
snphap.dir <- file.path(ROOT,if(args$cond) { "cond-snphap" } else { "snphap" },basename(args$d))
if(!file.exists(snphap.dir))
    dir.create(snphap.dir,recursive = TRUE)
  d <- genhaps(X,slist=hsnps,
               snps=DATA@snps,
               samples=df.samples,
               A1="A1",A2="A2",
               cols.copy=c("phenotype","country","PC1","PC2","PC3","PC4"),
               f.in=file.path(snphap.dir,"haps.in"),
               redo=TRUE)
 d$cc <- as.numeric(d$phenotype!="CONTROL")

#' haplotype analysis
library(mice)
RESULTS <- list()
for(ph in setdiff(DISEASES,c("RA.UK","RA.US"))) {
  covars <- if(ph %in% c("RA","COELIAC")) {
    use <- with(df.samples, phenotype %in% c("CONTROL",ph))
    covars <- c("country","PC1","PC2","PC3","PC4")
  } else {
    use <- with(df.samples, phenotype %in% c("CONTROL",ph) & country=="UK")
    covars <- c("PC1","PC2","PC3","PC4")
  }
  RESULTS[[ph]] <- model.mi(dir=snphap.dir,df.samples[use,], haps.pattern="haps.out2",covars=covars)
  RESULTS[[ph]]$disease <- ph
}

results <- lapply(RESULTS, function(x) {
  x <- x[grep("hap",rownames(x)),]
  x$hap <- factor(rownames(x),levels=rownames(x))
  return(x)
}) %>% do.call("rbind",.)
results$disease <- factor(results$disease)

#' haplotype frequencies
hfreq <- hapfreq(dir=snphap.dir,df.samples,haps.pattern="haps.out2",covars=c("country","phenotype"))
hfreq <- hfreq / matrix(rowSums(hfreq),nrow=nrow(hfreq),ncol=ncol(hfreq))
hfreq <- hfreq[!grepl("/NA",rownames(hfreq)), ]
write.csv(hfreq,file=myfile("haplotype-frequencies.csv"))

#' plot results
plots <- plotter(results,hsnps,hfreq)

library(ggbio)
p <- tracks(plots$odds,plots$freq + geom_hline(yintercept=1),plots$haps,heights=c(3,1,3),title=paste(REGION,NAME)) + xlim(0.75,nrow(plots$hlab)+1)
p
doplot(if(args$cond) { "cond-haplotypes" } else { "haplotypes" },p,width=10,type="pdf")

## #' diplotype analysis
## DRESULTS <- list()
## for(ph in setdiff(DISEASES,c("RA.UK","RA.US"))) {
##   if(ph=="RA") {
##     DRESULTS[[ph]] <-
##       diplotypes.mi(dir=snphap.dir,subset(df.samples,phenotype %in% c("CONTROL",ph)),
##                     haps.pattern="haps.out2",covars=c("country","PC1","PC2","PC3","PC4"))      
##   } else {
##     DRESULTS[[ph]] <- diplotypes.mi(dir=snphap.dir,subset(df.samples,phenotype %in% c("CONTROL",ph) & country=="UK"),
##                               haps.pattern="haps.out2",covars=c("PC1","PC2","PC3","PC4"))
##   }
##   DRESULTS[[ph]]$disease <- ph
## }

## dresults <- lapply(DRESULTS, function(x) {
##   x <- x[grep("dip",rownames(x)),]
##   x$dip <- factor(rownames(x),levels=rownames(x))
##   return(x)
## }) %>% do.call("rbind",.)
## dresults$disease <- factor(dresults$disease)

## dplots <- plotter(dresults,hsnps,hfreq,thr=0.03,y.min=-0.6,y.max=0.6)
## p <- tracks(dplots$odds,dplots$freq,dplots$haps,heights=c(3,1,4),title=paste(REGION,NAME)) + xlim(0.75,nrow(dplots$hlab)+1) 
## p

## doplot("diplotypes",p,width=10,height=10) 

