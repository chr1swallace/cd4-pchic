library(data.table)
library(magrittr)

## PMI
files <- list.files("/home/oliver/captureHIC_support/ICHIP/PMI_TAB_1",full=TRUE,
                    pattern="ATD|CEL|JIA|MS|PBC|RA|T1D")
PMI <- vector("list",length(files))
for(i in seq_along(files)) {
    pmi <- fread(files[[i]])
    setnames(pmi,names(pmi),c('chr','start','end','rs','r2','i.pos','maf','pval','ppi','delta'))
    PMI[[i]] <- pmi[,.(chr,start,end,rs,r2,pval,ppi)]
    PMI[[i]]$trait <- gsub(".*/|_IC.pmi","",files[[i]])
}
PMI <- do.call("rbind",PMI)

## NCVS
(load("~/FM/nsnps.RData"))
nsnps <- DATA
nsnps <- as.data.table(nsnps)
setkey(nsnps,region,trait)
nsnps[,post.nsnp:=sum(pp*nsnp),by=.(region,trait)]
nsnps <- nsnps[order(region,trait,pp,decreasing=TRUE),]
setkey(nsnps,region,trait)
nsnps <- unique(nsnps,by=c("region","trait"))

## Imputed manhattans
(load("~/FM/manhattan.RData"))

## MPPI
(load("~/FM/mppi.RData"))
MPPI <- as.data.table(MPPI)
setnames(MPPI,"var","snp")
MPPI$snp <- as.character(MPPI$snp)

setkey(nsnps,region,trait)
setkey(MPPI,region,trait)
MPPI <- merge(MPPI,nsnps[,.(trait,region,nsnp,post.nsnp)])

setkey(MANN,region,trait,snp)
setkey(MPPI,region,trait,snp)
MPPI <- merge(MPPI,MANN)

## merge
head(PMI)
head(MPPI)

## match traits
table(PMI$trait)
table(MPPI$trait)
MPPI$trait.orig <- MPPI$trait
MPPI[trait %in% c("COELIAC","ICOELIAC"), trait:="CEL"]
MPPI[trait=="GRAVES", trait:="ATD"]
MPPI[trait %in% c("RA","RA.UK","RA.US"), trait:="RA"]
table(MPPI$trait %in% PMI$trait)

## match snps
## table(MPPI$position %in% PMI$start)
## table(MPPI$position %in% (PMI$start+1))
## table(MPPI$position %in% PMI$end)
ss <- strsplit(MPPI$snp,"\\.")
rs <- sapply(ss,"[[",1)
pos <- sapply(ss,"[[",2) %>% as.numeric()
table(rs=rs %in% PMI$rs,pos=pos %in% PMI$start)
table(rs=rs %in% PMI$rs,pos=pos %in% PMI$end)
MPPI$rs <- rs

setkey(PMI,rs,trait)
setkey(MPPI,rs,trait)
m <- merge(MPPI,PMI)

setkey(m,region,trait)
m[,minpval:=min(pval), by=.(region,trait)]
m[,minp:=min(p), by=.(region,trait)]
m[,maxpval:=max(pval), by=.(region,trait)]
m[,maxp:=max(p), by=.(region,trait)]

with(m,table(ppi>0.5,mppi>0.5))
with(m,table(ppi>0.1,mppi>0.1))
with(m,table(ppi>0.01,mppi>0.01))
library(ggplot2)
ggplot(m,aes(x=ppi,y=mppi)) + geom_point() + geom_smooth() + geom_smooth(method="lm",se=FALSE) + facet_wrap(~trait)

dim(PMI)
dim(MPPI)
dim(m)

library(ggplot2)
library(ggbio)
ggplot(m,aes(x=ppi,y=mppi,col=as.factor(nsnp),lty=trait)) + geom_point() + facet_wrap(~region) + geom_smooth(method="lm",se=FALSE)

setkey(m,region,trait)
m2 <- m[minp<1e-6,]
m2[,cor:=cor(mppi,ppi),by=.(region,trait)]
m2 <- unique(m2)
 subset(m2,region=="10p-6030243-6169685")
ggplot(m2,aes(x=post.nsnp,y=cor,col=as.factor(trait))) + geom_point() + geom_hline(yintercept=0)

subset(m2,cor<0 & nsnp>1) ## only 2 regions: one is T1D/IL2RA, other T1D/INS
subset(m2,cor<0 & nsnp==1) ## 6 regions


