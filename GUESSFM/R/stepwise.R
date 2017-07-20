#!/usr/bin/Rscript
library(randomFunctions, quietly=TRUE)
library(annotSnpStats, quietly=TRUE)
library(magrittr)
source("DIRS.txt")

args <- getArgs(defaults=list(
                  file=file.path(ROOT,"FM-impute_qc/imputed-2q-204446380-204816382.RData"),
                  overwrite=FALSE,
                    #disease="JIA",
                  #region=75,
                  tag.r2=0.98, nsweep=50000, nsave=5000, nchains=5),
#                  tag.r2=0.6, nsweep=2000, nsave=2000, nchains=3),
                numeric=c("tag.r2","nsweep","nsave","nchains","region"))
# print(args)

# args$file <- "/stats/chrisw/FM-impute_qc/impute-14q-101290463-101328739.RData"

setwd("~/Projects/cd4chic/GUESSFM")

## files to read
d <- file.path(ROOT,"FM-impute_qc")
f.geno <- file.path(d,basename(args$file))
f.stepwise <- file.path(ROOT,"FM-stepwise",basename(args$file))
f.ss <- file.path(ROOT,"FM-ss/impdata",basename(args$file))
if(!file.exists(f.ss))
    stop("File not found: ",f.ss)

gd.root <- "/stats/chrisw/FM-GUESS"
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

## loop over phenotypes
phenotypes <- setdiff(unique(samples(DATA)$phenotype),c("CONTROL",NA))

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

if(!length(phenotypes))
  stop()

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

f.data <- function(ph) fd(file.path(gd,ph),"data.RData")

RESULTS <- structure(vector("list",length(phenotypes)),names=phenotypes)

for(ph in phenotypes) {
  message(ph)

  ## load data
(load(file=f.data(ph)))

  message(ph)
  message("Cases/controls: ", sum(cc.ph==1), " ", sum(cc.ph==0))
  message("SNPs: ",ncol(snp.data))

  X <- snp.data
  Y <- cc.ph
  
  ## matrix phenotypes
  if(!is.matrix(Y))
    Y <- matrix(Y,ncol=1)
  
  ## format snp data  
  p <- 0
  rs <- row.summary(X)
  
  use <- rs[,"Call.rate"]==1 & complete.cases(Y)
  if(!is.null(covars)) {
      use <- use & complete.cases(covars)
      p <- ncol(covars)
    }
  if(any(!use)) {
    use <- which(use)
    gX <- X[use,]
    gY <- Y[use]
    if(!is.null(covars))
      covars <- covars[use,]
  } else {
    gX <- X
    gY <- Y
  }
  n <- nrow(gX)
  cs <- col.summary(gX)
  use.cols <- which(!is.na(cs[,"z.HWE"]) & cs[,"MAF"]>0)
  m <- length(use.cols)
  if(m < ncol(gX))
    gX <- gX[,use.cols]

  if(!is.null(covars)) {
    if(!is.data.frame(covars) & !is.matrix(covars)) # deal with vectors, without messing up data.frames
      covars <- as.data.frame(as.matrix(covars))
    m0 <- glm(gY ~ ., data=covars,family="binomial")
    gY <- residuals(m0)
    family <- "gaussian"
  }
      
  message("using ",n," samples ",m," SNPs, ",p," covariates.")
  
  ## IF GUESS: what are the best regressors? (95,220?)
  mycond.best <- function(X,Y,best=NULL,stepwise.p.thr=1e-3,stepwise.max.predictors=NA, ...) {
 if(!is.na(stepwise.max.predictors) & length(best)>=stepwise.max.predictors)
   return(NULL)
 ## cs <- col.summary(X)
 ## X <- X[,cs[,"MAF"]>0.05]
 data <- data.frame(Y=Y, row.names=rownames(X))
 if(is.null(best)) {
   Xtest <- X
   cond <- single.snp.tests(phenotype=data$Y, snp.data=Xtest)
   p <- p.value(cond,1)
   stepwise.p.thr <- 1
  } else {
    data <- cbind(data,as(X[,best],"numeric"))
    LD <- ld(X,X[,best],stats="R.squared")
    maxLD <- apply(LD,1,max,na.rm=TRUE)
    drop <- unique(c(names(maxLD)[which(is.na(maxLD) | maxLD>0.5)],best))
    Xtest <- X[,setdiff(colnames(X),drop)]
    if(!ncol(Xtest))
      return(NULL)
    cond <- snp.rhs.tests(as.formula(paste("Y ~", paste(best, collapse="+"))),
                          snp.data=Xtest,
                          data=data,...) # binomial by default
    
    p <- p.value(cond)
  }
  pmin <- min(p, na.rm=TRUE)
  if(pmin<stepwise.p.thr) {
    newbest <- colnames(Xtest)[ which.min(p) ]
    cs <- col.summary(Xtest[,newbest])
    message(newbest,"\tMAF=",signif(cs[1,"MAF"],2),"\tp=",format.pval(pmin))
    return(list(snp=newbest,maf=cs[1,"MAF"],p=pmin))
  }
  return(NULL)
}

  best <- maf <- p <- NULL
  cs <- col.summary(X)
  wh <- which(!is.na(cs[,"z.HWE"]) & cs[,"MAF"]>0.005 & cs[,"Certain.calls"]>0.9)
  while(length(newbest <- mycond.best(X[use,wh], gY, best, family=family))) {
      maf <- c(maf,newbest[[2]])
      p <- c(p,newbest[[3]])
      best <- c(best,newbest[[1]])
  }

  RESULTS[[ph]] <- data.frame(trait=ph,snps=best,order=seq_along(best),maf=maf,p=p)
  
}

RESULTS <- do.call("rbind",RESULTS)

save(RESULTS,f.stepwise)
