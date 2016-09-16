#!/usr/bin/Rscript

library(randomFunctions)
library(magrittr)
args <- getArgs(defaults=list(file="./FM-impute/imputed-5p-489199-730271"))
library(snpStats)
library(annotSnpStats)

setwd(file.path(CD4CHIC.ROOT,"GUESSFM"))

message("!!! join-qc-imputed.R")
message(date())
message(getwd())
print(args)

## files to read
files <- list.files(dirname(args$file),pattern=paste0(basename(args$file),".gen.gz.gz$"))
if(length(files)==0)
    stop("no output file found for ", args$file)

## files to create
od <- file.path("./FM-impute_qc")
if(!file.exists(od))
  dir.create(od)
f.out <- file.path(od,paste0(basename(args$file),".RData"))
f.snps <- file.path(od,paste0(basename(args$file),".snps.gz"))
message("result will be saved to ",f.out)

d.data <- "./FM-impute"
d.pq <- file.path("./FM-aligned")

## start reading
source("R/impute_functions.R")
library(parallel)
if(interactive()) {
  options(mc.cores=20)
} else {
  options(mc.cores=1)
}

## sample info

chrarm <- gsub(".*/imputed-|-.*","",args$file)
chr <- sub("[pq]","",args$chrarm)
sample.data <- read.sinfo(d.pq,chrarm)

## impute results: info
f.info <- sub(".gz$","_info",files)
f.info.alt <- sub(".gen.gz.gz$",".nomono_info",files)
if(file.exists(file.path(d.data,f.info.alt)))
    f.info <- f.info.alt
message("reading impute SNP info from ",f.info)
impute.info <- read.imputeinfo(file.path(d.data,f.info))
rn <- make.unique(make.names(impute.info$rs_id))
rownames(impute.info) <- rn
message("SNPs found: ",nrow(impute.info))

## read samples
f.samples <- paste0(args$file,".sample") %>% sub("imputed","phased",.)
message("reading samples from ",f.samples)
samples <- read.table(f.samples,header=TRUE,as.is=TRUE)[-1,1]
message("samples found: ",length(samples))
sample.data <- sample.data[samples,]

## alleles
f.geno <- file.path(d.data,files)
f.geno.alt <- sub(".gen.gz.gz",".nomono_gen.gz",f.geno)
if(file.exists(f.geno.alt))
    f.geno <- f.geno.alt

message("reading alleles from ",f.geno)
alleles <- pipe(paste("zcat",f.geno,"| cut -d' ' -f 2,4,5")) %>%
  read.table(.,header=FALSE,sep=" ",col.names=c("rs_id","A1","A2"),as.is=TRUE)
head(alleles)
if(nrow(alleles) != nrow(impute.info))
  stop("SNP count mismatch between ",f.info," ",f.geno)
snp.info <- cbind(impute.info,alleles[,c("A1","A2")])

## genotypes
message("reading genotypes from ",f.geno)
X <- read.geno(f.geno,samples, alleles)

## check dimensions
message("read all genotypes\ndims: ")
dim(X)
message("sample dims:")
dim(sample.data)
message("snp dims:")
dim(snp.info)

if(nrow(sample.data)!=nrow(X))
    stop("sample count mismatch")
if(nrow(snp.info)!=ncol(X))
    stop("snp count mismatch")

## check sample order
rownames(sample.data) <- make.names(rownames(sample.data))
rownames(X) <- make.names(rownames(X))
if(!identical(rownames(sample.data),rownames(X)))
    stop("sample identifiers misordered")

## check SNPs
colnames(X) <- make.names(colnames(X))
rownames(snp.info) <- make.names(rownames(snp.info))
if(FALSE) { #debugging
  w <- which(colnames(X)!=rownames(snp.info))
  head(colnames(X)[w])
  head(snp.info[w,])
}
if(any(colnames(X)==".")) {
  w <- which(colnames(X)==".")
  snp.info[w,"rs_id"] <- paste0("chr",chrarm,":",snp.info[w,"position"])
  colnames(X)[w] <- rownames(snp.info)[w] <- snp.info[w,"rs_id"]
}

## start QC
csumm <- col.summary(X)
hwe.bycountry <- sapply(unique(sample.data$country), function(x)
                        col.summary(X[which(sample.data$affected==1 & sample.data$country==x),])[,"z.HWE"])
hwe.ctrl=apply(abs(hwe.bycountry),1,max,na.rm=TRUE)
print(ftable(rs.ok=snp.info$rs_id!=".", info.ok=snp.info$info>0.3, raf.nonzero=csumm[,"RAF"]>0, hwe.ctrl.ok=hwe.ctrl<5))
print(table(raf.nonzero=csumm[,"RAF"]>0, hwe.ctrl=hwe.ctrl<5))
  ## apply QC
use <- which(snp.info$rs_id!="." & snp.info$info>0.3 & csumm[ , "RAF"]>0 & hwe.ctrl<5)
DATA <- new("aSnpMatrix",
              .Data=X[,use],
            snps=snp.info[use,],
            samples=data.frame(samples=sample.data))
alleles(DATA) <- c("A1","A2")

save(DATA,file=f.out)
write.table(snps,file=f.snps)
