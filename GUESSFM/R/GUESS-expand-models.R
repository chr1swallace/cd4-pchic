#!/usr/local/Cluster-Apps/R/3.3.0/bin/Rscript
## follows on from GUESS-qc.R, may need some code from that inserted

library(randomFunctions)
library(annotSnpStats)
library(magrittr)
source("~/DIRS.txt")
setwd("~/Projects/cd4chic/GUESSFM")

                                        #library(GUESSFM)
library(devtools)
load_all("~/RP/GUESSFM")

## files
## 1p-2406887-2785671 y
## 15q-67414055-67469568 y
## 10q-59823894-60140709 y
## 13q-42834735-43101763 x
## 20p-1497197-1689461 x
## 21q-16698439-16843298 --
## 20p-1497197-1689461 y
## 21q-43810084-43887145 y
## 22q-21712414-22022005 y
## 2    204733868    204737096

## 2q-191873553-192007734
## calculating marginal SNP inclusion probabilities
## Error in .local(x, i, j, ..., drop) : 
##   No match for one or more column selections
## Calls: ld -> [ -> [ -> .local
## Execution halted

args <- getArgs(default=list(d=file.path(ROOT,"FM-GUESS/10p-6030243-6169685"),keep=FALSE))
args

source("R/GUESS-common.R")

library(ggplot2)
theme_set(theme_bw())
library(gridExtra)

dirs <- list.files(args$d,full=TRUE)
files <- unlist(lapply(dirs,list.files,full=TRUE,pattern="_features"))
                                        # files <- list.files(file.path(args$d,args$ph),full=TRUE
files <- grep("zzz|US|UK",files,invert=TRUE,value=TRUE)
files <- grep("ICOELIAC|T1D|RA|GRAVES",files,value=TRUE)
if(!length(files)) {
    message("No output files found")
    if(!interactive())
        q("no")
}
files <- sub("_sweeps_features.txt","",files)
                                        #files <- grep("_50000$",files,value=TRUE)
names(files) <- basename(dirname(files))

message("number of files found: ",length(files))
cat(files,sep="\n")

## ## for now
## if(length(wh <- which(names(files)=="ICOELIAC")))
##   files <- files[-wh]

## load guess 
dd <- lapply(files,read.snpmod)
names(dd) <- basename(dirname(files))

                                        # doplot("diffusion", plot_diffuse(dd))

## expand tags
TAGS <- lapply(files,function(f) {
    (load(file.path(dirname(f),"tags.RData")))
    return(tags)
})
names(TAGS) <- basename(dirname(files))


library(parallel)
NCORES <- if(interactive()) {1} else {8}
                                        #options(mc.cores=if(interactive()) {length(dd)} else {1})
options(mc.cores=NCORES) #if(system("hostname",intern=TRUE) %in% c("stats0","stats3")) {length(dd)} else {1})
                                        #system(sprintf('taskset -p 0xffffffff %d', Sys.getpid()))

## trim least probable models
## trim <- function(x,thr=1e-5) {
##     use <- x@models$PP>thr
##     x@models <- x@models[use,]
##     x@model.snps <- x@model.snps[use]
##     x
## }
## dd <- mclapply(dd,trim)

dx <- mcmapply(expand.tags,dd,TAGS,SIMPLIFY=FALSE)
                                        # best.models(dd)
bm <- best.models(dx)

##' FIRST, STOP HERE IF NULL MODEL HAS > 50% POSTERIOR ACROSS ALL TRAITS
dropmodels <- function(x,minpp=0.5) {
    if(!any(x$str==""))
        return(FALSE)
    w <- which(x$str=="")
    x[w,"PP"]
}

nullpp <- lapply(bm, dropmodels) %>% unlist()
message("null model posterior")
nullpp
message("null model >50% of posterior?")
print(drop <- nullpp>0.5)

if(all(drop)) {
    system(paste0("touch ",file.path(args$d,"skip")))
    if(!interactive())
        q("no")
}

## GIVEN SOME TRAITS HAVE MODELS OF INTEREST, DROP ONLY THOSE THAT HAVE REALLY NO SUPPORT
drop <- nullpp>0.8

if(any(drop)) {
    drop <- names(drop)[drop]    
    message("dropping ",paste(drop,collapse=" "))
    dd <- dd[setdiff(names(dd),drop)]
    dx <- dx[setdiff(names(dd),drop)]
    dirs <- dirs[-grep(paste(drop,collapse="|"),dirs)]
    files <- files[-grep(paste(drop,collapse="|"),files)]
}

## best.snps(dd)
## best.snps(dx)


## TODO: LOOK AT BEST SNPS - WHAT IS THEIR QC?

## v <- unique(unlist(lapply(best.snps(dd),rownames)))
## col.summary(DATA[,v])
## p <- samples(DATA)
## col.summary(DATA[which(p$phenotype=="CONTROL"),v])
                                        #for(p in unique(p$phenotype)) {
#if(FALSE) {
    
    ## load DATA
dfiles <- unlist(sapply(dirs,list.files,pattern="data.RData",full=TRUE))
names(dfiles) <- basename(names(dfiles))
dfiles <- dfiles[names(files)]
    DATA <- lapply(dfiles, function(f) { (load(f)); return(as(snp.data[cc.ph==0,],"SnpMatrix"))})
    
    SP <- mcmapply(snp.picker,dd,DATA)
    SPx <- mcmapply(snp.picker,dx,DATA)
    
    ## pdf(file.path(plot.dir,"snppicker.pdf"),height=8,width=8)
    ## for(i in seq_along(SP))
    ##   grid.newpage()
    ##   p <- plot(SP[[i]],do.plot=FALSE)
    ##   print(p + ggtitle(names(files)[[i]]))
    ## dev.off()
    
    
    ## groups <- as(SP[[1]],"groups")
    ## for(i in 2:length(SP))
    ##   groups <- union(groups,as(SP[[i]],"groups"))
    
    library(reshape)
    library(plyr)
    
    SPx
    groupsx <- as(SPx[[1]],"groups")
    if(length(SPx)>1)
        for(i in 2:length(SPx))
            groupsx <- union(groupsx,as(SPx[[i]],"groups"))
        
    ## pattern.plot(dd,groups)
    doplot("pattern-approx", pattern.plot(dx,groupsx))
#}

## refits
f.geno <- file.path(ROOT,"FM-impute_qc",paste0("imputed-",basename(args$d),".RData"))
(load(f.geno))
colnames(samples(DATA)) <- sub("samples.","",colnames(samples(DATA)))
message("region loaded from ",f.geno)
message(ncol(DATA)," SNPs in region")
phenotype(DATA) <- "affected"

## manhattans
file.manhattan <- file.path(args$d,"manhattan.RData")
## if(file.exists(file.manhattan)) {
##     (load(file.manhattan))
##     if("position" %in% names(MANHATTAN[[1]]))
##         next
## }

MANHATTAN <- vector("list",length(files))
names(MANHATTAN) <- names(files)
for(ph in names(files)) {
    message("single.snp.tests for ",ph)
    dfile <- file.path(dirname(files[[ph]]),"data.RData")
    tfile <- file.path(dirname(files[[ph]]),"tags.RData")
    (load(dfile))
    (load(tfile))
    snp.data <- snp.data[,snps(tags)] # only those SNPs which passed QC
    tmp <- single.snp.tests(phenotype=cc.ph,snp.data=snp.data,stratum=strat)
    MANHATTAN[[ph]] <- data.frame(snp=tmp@snp.names,p=p.value(tmp,1),trait=ph)
    ## if(ph=="RA") {
    ##     for(s in c("UK","US")) {
    ##         phs <- paste(ph,s,sep=".")
    ##         message("single.snp.tests for ",phs)
    ##         tmp <- single.snp.tests(phenotype=cc.ph,snp.data=snp.data,
    ##                                 subset=strat==s)
    ##         MANHATTAN[[phs]] <- data.frame(snp=tmp@snp.names,p=p.value(tmp,1),trait=phs)
    ##     }
    ## } 
}

## proxies <- c("rs231779",
## "rs231775",
## "rs231723",
## "rs231724",
## "rs11571315",
## "rs231770",
## "rs231764",
## "rs926169",
## "rs1024161",
## "rs231763",
## "rs11571292",
## "rs4675372")

## kk2 <- kk[ grep(paste(proxies,collapse="|"),kk$snp),]
## kk2[order(kk2$position),c("snp_id","rs_id","position","info","p")]

## add snp info
ss <- DATA@snps[,c("snp_id","rs_id","position","exp_freq_a1","info","certainty","A1","A2")]
ss$snp <- make.names(ss$rs_id)

for(ph in names(MANHATTAN))
    MANHATTAN[[ph]] %<>% merge(.,ss,by="snp",all.x=TRUE)

save(MANHATTAN,file=file.manhattan)
## }


library(speedglm)
library(parallel)
options(mc.cores=NCORES) #if(system("hostname",intern=TRUE) %in% c("stats0","stats1")) { 20 } else { 1 })

tmp <- best.models(dx,cpp.thr=0.9)
if(length(tmp)>1) {
  union <- model.union(tmp[[1]]$str,tmp[[2]]$str,detail=FALSE)
  if(length(dx)>2)
    for(j in 3:length(dx))
      union <- model.union(union,tmp[[j]]$str)
  length(union)
} else {
  union <- tmp[[1]]$str
}

fits <- structure(vector("list",length(files)), names=names(files))

N <- vector("list",length(files))
names(N) <- names(files)
for(ph in names(files)) {
  dfile <- file.path(dirname(files[[ph]]),"data.RData")
  (load(dfile))  
  snp.data <- snp.data[,TAGS[[ph]]@.Data]
  fits[[ph]] <- abf.calc(y=cc.ph,x=snp.data,models=union,family="binomial",
                             snp.data=as(DATA,"SnpMatrix"),
                             q=covars)
  ## if(ph=="RA") {
  ##   for(s in c("UK","US")) {
  ##     wh <- which(strat==s)
  ##     phs <- paste(ph,s,sep=".")
  ##     fits[[phs]] <- abf.calc(y=cc.ph[wh],x=snp.data[wh,],models=union,family="binomial",
  ##                                 snp.data=as(DATA,"SnpMatrix"),
  ##                                 q=covars[wh,])
  ##   }
  ## } 
}

allph <- names(files)
## if("RA" %in% names(dx)) {
##   ph <- "RA"
##   dfile <- file.path(dirname(files[[ph]]),"data.RData")
##   (load(dfile))
##   allph <- c(allph,paste("RA",c("UK","US"),sep="."))
## } 


## make snpmods
SM2 <- lapply(fits,function(x) abf2snpmod(x,expected=2,overdispersion=1))
save(SM2,file=file.path(args$d,"snpmod.RData"))

## LD <- ld(DATA,DATA[,"rs231779.204734487.C.T",drop=FALSE],stat="R.squared")
## LD <- data.frame(snp=rownames(LD),r2=LD[,1])
## x <- merge(x,LD,by="snp",all.x=TRUE)
## ggplot(x, aes(x=position,y=-log10(p),col=1-r2)) + geom_point() + facet_grid(trait ~ .)
## ggsave("ctla4-manhattans.pdf",height=8,width=10)
## paste0("scp stats0:",getwd(),"/ctla4-manhattans.pdf"," .")

SM2
best.snps(SM2)
best.models(SM2)

nsnp <- pp.nsnp(SM2,expected=2)
doplot("nsnps", nsnp@plot, display=TRUE)

v <- setdiff(unique(unlist(lapply(best.snps(SM2),rownames))),"1")
if(length(v))
    print(col.summary(DATA[,v]))
p <- samples(DATA)

## if(args$name=="IL2RA") {

## v <- setdiff(unique(unlist(lapply(best.snps(SM2),rownames))),"1") 
## col.summary(DATA[,v])
## p <- samples(DATA)
## col.summary(DATA[which(p$phenotype=="CONTROL" & p$country=="UK"),v])
## col.summary(DATA[which(p$phenotype=="COELIAC" & p$country=="UK"),v])
## col.summary(DATA[which(p$phenotype=="RA" & p$country=="UK"),v])
## col.summary(DATA[which(p$phenotype=="RA" & p$country=="US"),v])
## col.summary(DATA[which(p$phenotype=="CONTROL" & p$country=="UK"),v])
## col.summary(DATA[which(p$phenotype=="CONTROL" & p$country=="US"),v])

## v <- row.names(best.snps(SM2$RA.UK,0.05))
## p <- samples(DATA)
## col.summary(DATA[which(p$phenotype=="CONTROL" & p$country=="UK"),v])
## col.summary(DATA[which(p$phenotype=="COELIAC" & p$country=="UK"),v])
## col.summary(DATA[which(p$phenotype=="RA" & p$country=="UK"),v])
## col.summary(DATA[which(p$phenotype=="RA" & p$country=="US"),v])
## col.summary(DATA[which(p$phenotype=="CONTROL" & p$country=="US"),v])

## }


## if(args$name=="CTLA4") {
##   for(ph in setdiff(unique(p$phenotype),NA)) {
##     message(ph)
##     print(col.summary(DATA[which(p$phenotype==ph & p$country=="UK"),v])[,c("Calls","Certain.calls","RAF","z.HWE")])
##   }

##   ## rs74864202.204539270.C.A bad in JIA
## }


## group snps
dfiles <- unlist(sapply(dirs,list.files,pattern="data.RData",full=TRUE))
names(dfiles) <- basename(dirname(dfiles))
DATA <- lapply(dfiles, function(f) { (load(f)); return(as(snp.data[cc.ph==0,],"SnpMatrix"))})
## if("RA" %in% allph) {
##   (load(dfiles[["RA"]]))
##   nm <- grep("^RA.",names(fits),value=TRUE)
##   for(n in nm) {
##     DATA[[n]] <- as(snp.data[cc.ph==0 & strat==sub("RA.","",n), ], "SnpMatrix")
##   }
## }
#DATA[["RA2"]] <- DATA[["RA"]]

## dropsnps <- function(d,snps) {
##     d@snps <- subset(d@snps,!(var %in% snps))
##     moddrop <- grep(paste(snps,collapse="|"),d@model.snps)
##     if(length(moddrop)) {
##         d@model.snps <- d@model.snps[-moddrop]
##         d@models <- d@models[-moddrop,]
##     }
##     return(d)
## }

## if(args$d=="/stats/chrisw/FM-GUESS/22q-21712414-22022005") {
##     SM2[["JIA"]] <- dropsnps (SM2[["JIA"]], "rs76690447.21998012.C.A")
## }
## if(args$d=="/stats/chrisw/FM-GUESS/1p-2406887-2785671") {
    
## }
## a <- d@snps
##   a$var <- as.character(a$var)
##   a <- subset(a,var!="1") # ignore null model
## setdiff(a$var,colnames(data))

if(length(DATA)>1 & "JIA" %in% names(DATA)) {
    ok <- c("RA.UK","T1D","MS","GRAVES","COELIAC")
    if(length(wh <- which(ok %in% names(DATA))))
        DATA[["JIA"]] <- DATA[[ ok[wh[1]] ]]
}

## pattern plots
SP2 <- mapply(snp.picker,SM2,DATA)
save(SP2,file=file.path(args$d,"snp-picker.RData"))
groups2 <- as(SP2[[1]],"groups")
if(length(SP2)>1)
  for(i in 2:length(SP2))
    groups2 <- GUESSFM::union(groups2,as(SP2[[i]],"groups"))

doplot("pattern-refit",pattern.plot(SM2,groups2))

message("--END--")


# CT60 = rs3087243
# +49 G/A = rs231775

## > best.snps(SM2)
## $COELIAC
##                         Marg_Prob_Incl                     var
## rs3087243.204738919.G.A      0.4269163 rs3087243.204738919.G.A
##                                        rownames
## rs3087243.204738919.G.A rs3087243.204738919.G.A

## $GRAVES
##                          Marg_Prob_Incl                      var
## rs3087243.204738919.G.A       0.2633483  rs3087243.204738919.G.A
## rs11571297.204745003.T.C      0.2098508 rs11571297.204745003.T.C
##                                          rownames
## rs3087243.204738919.G.A   rs3087243.204738919.G.A
## rs11571297.204745003.T.C rs11571297.204745003.T.C

## $RA
##                          Marg_Prob_Incl                      var
## rs62184035.204661674.A.G      0.7465553 rs62184035.204661674.A.G
## rs3087243.204738919.G.A       0.2215348  rs3087243.204738919.G.A
## rs1427676.204741166.C.T       0.1960494  rs1427676.204741166.C.T
## rs231729.204743801.T.A        0.1931947   rs231729.204743801.T.A
## rs231725.204740675.G.A        0.1094215   rs231725.204740675.G.A
## rs231726.204740866.C.T        0.1016959   rs231726.204740866.C.T
##                                          rownames
## rs62184035.204661674.A.G rs62184035.204661674.A.G
## rs3087243.204738919.G.A   rs3087243.204738919.G.A
## rs1427676.204741166.C.T   rs1427676.204741166.C.T
## rs231729.204743801.T.A     rs231729.204743801.T.A
## rs231725.204740675.G.A     rs231725.204740675.G.A
## rs231726.204740866.C.T     rs231726.204740866.C.T

## $T1D
##                         Marg_Prob_Incl                     var
## rs4673266.204634569.A.T      0.6411179 rs4673266.204634569.A.T
## rs231775.204732714.A.G       0.3319982  rs231775.204732714.A.G
## rs231724.204739823.A.G       0.2841487  rs231724.204739823.A.G
## rs231723.204739781.A.G       0.2695615  rs231723.204739781.A.G
##                                        rownames
## rs4673266.204634569.A.T rs4673266.204634569.A.T
## rs231775.204732714.A.G   rs231775.204732714.A.G
## rs231724.204739823.A.G   rs231724.204739823.A.G
## rs231723.204739781.A.G   rs231723.204739781.A.G

################################################################################

################################################################################

#' copy all output back to main project directory
## if(interactive())
##    system(paste("rsync -av ",plot.dir," chrisw@salthouse:~/Projects/FineMapping/output"))
