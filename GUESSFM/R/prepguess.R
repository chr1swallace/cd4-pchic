#!/usr/local/Cluster-Apps/R/3.3.0/bin/Rscript
library(randomFunctions, quietly=TRUE)
library(annotSnpStats, quietly=TRUE)
library(magrittr)
source("~/DIRS.txt")

args <- getArgs(defaults=list(
                  file=file.path(ROOT,"FM-impute_qc/imputed-12q-68429119-68510072.RData"),
                  overwrite=FALSE,
                    #disease="JIA",
                  #region=75,
                  tag.r2=0.98, nsweep=50000, nsave=5000, nchains=5),
#                  tag.r2=0.6, nsweep=2000, nsave=2000, nchains=3),
                numeric=c("tag.r2","nsweep","nsave","nchains","region"))
# print(args)

# args$file <- "/stats/chrisw/FM-impute_qc/impute-14q-101290463-101328739.RData"

setwd("~/Projects/cd4chic/GUESSFM")


## if(args$file=="/stats/chrisw/FM-impute_qc/impute-2q-204446380-204816382.RData" & !interactive())
##   stop("not CTLA4!")


## files to read
d <- file.path(ROOT,"FM-impute_qc")
## (load(file.path(d,"common","iChipFineMappingRegionsB37.RData")))
## as.data.frame(ichip.regions.37)

f.geno <- file.path(d,basename(args$file))
f.ss <- file.path(ROOT,"FM-ss/impdata",basename(args$file))
if(!file.exists(f.ss))
    stop("File not found: ",f.ss)

## create vector, snps.drop
source("R/drop-imputed-snps.R")

## files to create
##  library(devtools, quietly=TRUE)
## install_github("chr1swallace/GUESSFM")
library(GUESSFM)
#load_all("~/RP/GUESSFM",quiet=TRUE)
gd.root <- file.path(ROOT,"FM-GUESS")
file.root <- gsub("impute-|imputed-|.imputed.RData|.RData","",basename(f.geno))

gd <- file.path(gd.root,file.root)
if(!file.exists(gd))
  dir.create(gd,recursive=TRUE)
fd <- function(d,f) {
  if(!file.exists(d))
    dir.create(d,recursive=TRUE)
  file.path(d,f)
}

## load data and phenotypes
message(f.geno)
print(load(f.geno))
colnames(samples(DATA)) <- sub("samples.","",colnames(samples(DATA)))

message("region loaded from ",f.geno)
message(ncol(DATA)," SNPs in region")
phenotype(DATA) <- "affected"

## drop bad snps
wh <- which(annotSnpStats::snps(DATA)$rs_id %in% snps.drop |
            annotSnpStats::snps(DATA)$snp_id %in% snps.drop |
            rownames(annotSnpStats::snps(DATA)) %in% snps.drop)
message("dropping ",length(wh)," SNPs listed in snps.drop")
if(length(wh)) {
  DATA <- DATA[,-wh]
}

snps.drop3 <- rownames(DATA@snps)[(DATA@snps$info<0.3 & DATA@snps$type==0) |
                                  (DATA@snps$certainty<0.98 & DATA@snps$type==0) ]
wh <- which(annotSnpStats::snps(DATA)$rs_id %in% snps.drop3 |
            annotSnpStats::snps(DATA)$snp_id %in% snps.drop3 |
            rownames(annotSnpStats::snps(DATA)) %in% snps.drop3)
message("dropping ",length(wh)," SNPs with info<0.3 or certainty<0.98")
if(length(wh)) {
  DATA <- DATA[,-wh]
}

## non-missing phenotype
cc <- samples(DATA)$affected - 1
use <- complete.cases(cc)
cc <- cc[which(use)]
DATA <- DATA[which(use),]
message("subsetting to ",length(which(use))," SAMPLEs with phenotype information")
print(table(cc))

with(samples(DATA),table(file,phenotype))

## loop over phenotypes
phenotypes <- setdiff(unique(samples(DATA)$phenotype),"CONTROL")
if(!length(phenotypes)){
  paste("touch",file.path(gd,"skip")) %>% system()
}

library(fastcluster,quietly=TRUE) # make hclust faster: O(n^3) -> O(n^2)
library(parallel,quietly=TRUE)

## what is min p?  Only go ahead if mip < 10-5
(load(f.ss))

transph <- function(ph) { # go from Nick's names to mine
  ph <- tolower(ph)
  ph <- sub("coeliac","celiac",ph)
  ph <- sub("graves","atd",ph)
  return(ph)
}

pvars <- paste("p",transph(phenotypes),sep=".")
print(minp <- apply(results[,pvars],2,min,na.rm=TRUE))
phenotypes <- phenotypes[ minp < 1e-5 ]
sapply(phenotypes, function(ph) message("INTENDED ",basename(args$file)," ",ph))

if(!length(phenotypes)) {
  message()
  message("--------------------------------------------")
  message("no significant assoc to the region, quitting")
  message("--------------------------------------------")
  message()
  q("no")
}

## quick check
tagfiles <- file.path(d,phenotypes,"tags.RData")
if(!interactive() && all(sapply(tagfiles,file.exists))) {
  message()
  message("----------------------------------------")
  message("all output files already exist, quitting")
  message("----------------------------------------")
  message()
  q("no")
}


## look at differential certain call rates
PH <- c("CONTROL",phenotypes)
CS <- lapply(PH,function(ph) {
  col.summary(DATA[ which(samples(DATA)$phenotype==ph), ])
})
todrop <- function(x,y) {
  abs(log(x)-log(y))>0.05 | (y>0.99 & x<0.98) | (x>0.99 & y<0.98)
}
w1 <- which(PH=="CONTROL")
keep <- rep(TRUE,ncol(DATA))
for(ph in setdiff(PH,"CONTROL")) {
  w2 <- which(PH==ph)
  x <- CS[[w1]][,"Certain.calls"]
  y <- CS[[w2]][,"Certain.calls"]
  keep <- keep & !todrop(x,y)  
}
if(any(!keep)) {
  message("dropping ",sum(!keep)," SNPs with differential certain calls")
  DATA <- DATA[,which(keep)]
}

## regions/phenotypes to drop
rp <- read.table("drop-regions-phenotypes.txt",as.is=TRUE,header=TRUE)
wh <- which(rp$region==basename(gd) & rp$ph %in% phenotypes)
if(length(wh)) {
  message("dropping region/phenotypes from drop-regions-phenotypes.txt: ",length(wh))
  print(rp[wh,])
  phenotypes <- setdiff(phenotypes,rp[wh,"ph"])
}
if(!length(phenotypes))
  stop()

## add pcs
f.samples <- file.path(ROOT,"FM-pca/PCS+.RData")
(load(f.samples))
rownames(pcs) <- make.names(rownames(pcs))
sinfo <- samples(DATA)
sinfo <- pcs[rownames(sinfo),]

message("Running GUESS for ",length(phenotypes)," diseases: ",paste(phenotypes,collapse=", "))
  options(scipen=1000000000)

coms <- structure(rep("",length(phenotypes)),names=phenotypes)
for(ph in phenotypes) {
  message(ph)

  f.tags <- function(ph) fd(file.path(gd,ph),"tags.RData")
  f.data <- function(ph) fd(file.path(gd,ph),"data.RData")
  f.par <- function(ph) fd(file.path(gd,ph), "par.xml")
  if(file.exists(f.par(ph)) & !args$overwrite) {
    message("output file already exists, skipping: ",f.par(ph))
    next
  }

if(ph %in% c("RA","ICOELIAC")) {
  wh.DATA <- which(sinfo$phenotype %in% c(ph,"CONTROL") # &
  #                   sinfo$outlier.code==0
)
} else {
  wh.DATA <- which(sinfo$phenotype %in% c(ph,"CONTROL") &
   #                sinfo$outlier.code==0 &
                   sinfo$country=="UK")
}
cc.ph <- cc[wh.DATA]
snp.data <- as(DATA,"SnpMatrix")[wh.DATA,]

  ## PCs
  covars <- sinfo[wh.DATA,c("PC1","PC2","PC3","PC4")]
  fix <- which(sapply(covars,class)=="matrix")
  if(length(fix)) {
    for(j in fix) 
      covars[,j] <- as.vector(covars[,j])
  }
  
## stratify if RA/ICOELIAC
strat <- factor(sinfo[wh.DATA,"country"])
if(length(levels(strat))>1)
  covars <- cbind(covars,strat=strat)

  ## save data
save(snp.data,cc.ph,strat,covars,file=f.data(ph))

## tags
  if(args$tag.r2<1) {
    message("Preparing to tag.")
    message("Input matrix has ",ncol(snp.data)," SNPs.")
    cs <- col.summary(snp.data)
    cs0 <- col.summary(snp.data[cc.ph==0,])
    cs1 <- col.summary(snp.data[cc.ph==1,])
    wh <- which(is.na(cs0[,"z.HWE"]) | is.na(cs1[,"z.HWE"]) |
                cs0[,"MAF"]<0.01 |
                cs0[,"Call.rate"]<0.99 | cs1[,"Call.rate"]<0.99 |
                cs0[,"Certain.calls"]<0.75 | cs1[,"Certain.calls"]<0.75 |
                abs(cs0[,"z.HWE"])>4)
    if(length(wh)) {
      message("Dropping ",length(wh)," SNPs with |z.HWE|>5, MAF < 0.01 in controls or call rate <0.99 in cases or controls")
      snp.data <- snp.data[,-wh]
    }
    tags <- tag(snp.data, method="single", tag.threshold=args$tag.r2)
    snp.data <- snp.data[, unique(tags(tags))]
    message("after tagging at ",args$tag.r2,", matrix now has ",ncol(snp.data)," SNPs.")
    save(tags, file=f.tags(ph))
  }
  
  message(ph)
  message("Cases/controls: ", sum(cc.ph==1), " ", sum(cc.ph==0))
  message("SNPs: ",ncol(snp.data))

  ## command=paste("qCom.sh -L", file.path(gd,ph,"2log"), "/home/chrisw/local/bin/GUESS")
  ## command=paste("/home/chrisw/local/bin/GUESS")
  coms[[ph]] <- run.bvs(X=snp.data,Y=cc.ph,covars=covars,
                       gdir=file.path(gd,ph),nsweep=args$nsweep, family="binomial",tag.r2=NA,
                       nsave=args$nsave,nchains=args$nchains,nexp=2,run=FALSE)

  ## if(ph=="RA") { # do RAUK alone
  ##   (load(f.data(ph)))
  ##   use <- strat=="UK"
  ##   snp.data <- snp.data[use,]
  ##   cc.ph <- cc.ph[use]
  ##   covars <- covars[,1:(ncol(covars)-1)]
    
    
  ## }
  
}

  cat(coms,file=tempfile("runguess",tmpdir=".",fileext="sh"),sep="\n")
