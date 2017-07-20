library(GEOquery)
library(ggplot2)
source("common.R")

f <- "GSE60680.RData"
if(!file.exists(f)) {
  message("retrieving GSE60680 from GEO")
  gsm <- getGEO("GSE60680")
  save(gsm,file=f)
} else {
  message("loading local copy of GSE60680")
  load(f)
}

pd <- lapply(gsm,pData)
lapply(pd,colnames)
lapply(pd, function(x) head(x[,grep("charact",colnames(x))],20))
fd <- lapply(gsm,fData)
lapply(fd,colnames)
lapply(fd,function(x) x[100,])
lapply(pd, function(x) head(x[,grep("charact",colnames(x))],20))

f1 <- fData(gsm[[2]])
f2 <- fData(gsm[[5]])
p1 <- pd[[2]][,c("title","source_name_ch1",grep("charact",colnames(pd[[2]]),value=TRUE))]
p2 <- pd[[5]][,c("title","source_name_ch1",grep("charact",colnames(pd[[2]]),value=TRUE))]

library(preprocessCore)
use <- f1$CONTROL_TYPE=="FALSE"
x1 <- log2(exprs(gsm[[2]]))[use,]
f1 <- f1[use,]

use <- f2$CONTROL_TYPE=="FALSE"
x2 <- log2(exprs(gsm[[5]]))[use,]
f2 <- f2[use,]

## FOLD CHANGE
rn <- intersect(rownames(x1),rownames(x2))
#x <- cbind(x1[rn,],x2[rn,])
x <- x1[rn,]
x <- normalize.quantiles(x)
p <- p1 #rbind(p1,p2)
cn <- intersect(colnames(f1),colnames(f2))
f <- f1[rn,] #rbind(f1[rn,cn],f2[rn,cn])

library(limma)
library(data.table)
library(magrittr)
p <- as.data.table(p)
setnames(p,c("title","skip","cell","time"))
p[,skip:=NULL]
ss <- strsplit(as.character(p$title), "-")
p[,sample.id:=sapply(ss,"[[",1)]
p[,time:=gsub("time: ","",time) %>% factor(.,levels=c("0h","6h","24h","3d","6d","8d"))]
p[,cell:=gsub("cell type: ","",cell)]
p[,variable:=colnames(x1)]
stim.times <- setdiff(unique(p$time), "0h")
stim.cells <- setdiff(unique(p$cell), "Naive T")
p$Th2 <- p$cell=="Th2"


mods <- fread(file.path(CD4CHIC.EXTDATA,"gene-modules.csv"))
#mods <- as.data.frame(modules)
mods <- mods[!duplicated(mods$geneSymbol),]
#rownames(mods) <- mods$geneSymbol
m <- match(as.character(f$GENE_SYMBOL),mods$geneSymbol)
f$module <- NA
f[!is.na(m),"module"] <- mods[ m[!is.na(m)], ]$all_12

Th1 <- which(p$cell %in% c("Th1","Naive T"))
design <- with(p[Th1,], model.matrix(~ sample.id + time))
efit.Th1 <- lmFit(x[,Th1],design) %>% eBayes()
efit.Th1$genes <- f[,c("ID","GENE_SYMBOL","module"),drop=FALSE]
Th2 <- which(p$cell %in% c("Th2","Naive T"))
design <- with(p[Th2,], model.matrix(~ sample.id + time))
efit.Th2 <- lmFit(x[,Th2],design) %>% eBayes()
efit.Th2$genes <- f[,c("ID","GENE_SYMBOL","module"),drop=FALSE]

RESULTS <- list(Th1.6h=topTable(efit.Th1,coef="time6h",p.value=1,number=Inf),
                Th1.24h=topTable(efit.Th1,coef="time24h",p.value=1,number=Inf),
                Th1.3d=topTable(efit.Th1,coef="time3d",p.value=1,number=Inf),
                Th1.6d=topTable(efit.Th1,coef="time6d",p.value=1,number=Inf),
                Th1.8d=topTable(efit.Th1,coef="time8d",p.value=1,number=Inf),
                Th2.6h=topTable(efit.Th2,coef="time6h",p.value=1,number=Inf),
                Th2.24h=topTable(efit.Th2,coef="time24h",p.value=1,number=Inf),
                Th2.3d=topTable(efit.Th2,coef="time3d",p.value=1,number=Inf),
                Th2.6d=topTable(efit.Th2,coef="time6d",p.value=1,number=Inf),
                Th2.8d=topTable(efit.Th2,coef="time8d",p.value=1,number=Inf))
lapply(RESULTS,function(x) table(x$module))
for(nm in names(RESULTS))
    RESULTS[[nm]]$test <- nm

RESULTS <- do.call("rbind",RESULTS)
table(RESULTS$module)

RESULTS <- as.data.table(RESULTS)
RESULTS[,avg.logFC:=median(logFC),by=c("module","test")]
for(i in list(0.05,0.1,0.25,0.5,0.75,0.9,0.95)) {
    nm=paste0("q",round(100*i))
    RESULTS[,c(nm) := list(quantile(logFC,i)),
            by=c("module","test")]
}

RESULTS <- unique(RESULTS,by=c("module","test"))
RESULTS[,cell:=sub("\\..*","",test)]
RESULTS[,time:=factor(sub(".*\\.","",test), levels=c("0h","6h","24h","3d","6d","8d"))]
t0 <- RESULTS[time=="6h",]
t0[,avg.logFC:=0]
for(i in list(0.05,0.1,0.25,0.5,0.75,0.9,0.95)) {
    nm=paste0("q",round(100*i))
    t0[,c(nm) := list(0)]
}
t0[,time:="0h"]
RESULTS <- rbind(t0,RESULTS)

colScale <- make.colScale(unique(RESULTS$module))
ggplot(RESULTS,aes(x=time,y=avg.logFC,col=module,group=module)) + geom_point() + geom_path() +
  geom_path(aes(y=q25),linetype="dotted") + geom_path(aes(y=q75),linetype="dotted") +
facet_grid(cell ~ .) + colScale

## our results
short <- fread(file.path(CD4CHIC.EXTDATA,"microarray-diffexpr.csv"))
mods <- as.data.table(mods)
short <- merge(short,mods[,.(id,all_12)],by.x="ens.gene.id",by.y="id")
setnames(short,"all_12","module")
short[,test:=paste0("total.",time,"h")]
short[,time:=paste0(time,"h")]
short[,avg.logFC:=median(logFC),by=c("module","test")]
for(i in list(0.05,0.1,0.25,0.5,0.75,0.9,0.95)) {
    nm=paste0("q",round(100*i))
    short[,c(nm) := list(quantile(logFC,i)),
            by=c("module","test")]
}
short <- unique(short,by=c("module","test"))
short[,cell:="short"]
short0 <- short[time=="2h",]
short0[,time:="0h"]
short0[,avg.logFC:=0]
short <- rbind(short0,short)

kk <- rbind(RESULTS[time!="0h",.(module,cell,time,avg.logFC,q25,q75)],
            short[,.(module,cell,time,avg.logFC,q25,q75)])
kk <- kk[!is.na(module) & cell!="Th1",]
kk[time=="0h",q25:=0]
kk[time=="0h",q75:=0]
days <- grepl("d",kk$time)
kk$timen <- as.numeric(sub("[hd]","",kk$time)) + ifelse(days,24,0)
kk$time <- factor(as.character(kk$time),levels=c("0h","2h","4h","6h","21h","24h","3d","6d","8d"))
kk <- kk[order(time,cell),]
kk$module <- as.factor(kk$module)
colScale <- make.colScale(unique(kk$module))
w <- 0.04
nmod <- length(unique(kk$module))
kk$cell <- ifelse(kk$cell=="short",
                  "Total CD4, short timecourse",
                  paste("Naive CD4, long timecourse"))
kk$cell <- factor(kk$cell,levels=rev(sort(unique(kk$cell))))
ggplot(kk,
       aes(##x=timen,
         ##x=ifelse(time=="0h",1,as.numeric(time)+w*as.numeric(module)-nmod*w/2),
         x=as.numeric(time) + as.numeric(time=="6h") * ifelse(cell=="Naive CD4, long timecourse",0.05,-0.05),
           y=avg.logFC,col=module,pch=cell)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) + geom_path(aes(linetype=cell)) +
  geom_errorbar(aes(ymin=q25,ymax=q75),alpha=0.5,width=.2) +
  ##geom_path(aes(y=q25),linetype="dotted") + geom_path(aes(y=q75),linetype="dotted") +
  ##  facet_grid(cell ~ .) +
  facet_grid(module~.,scales="free") +
  colScale +
  ylab("Average log FC") +
  scale_shape_manual(values=c(15,5)) +
  scale_x_continuous("Time",
                     breaks=sort(unique(as.numeric(kk$time))),
                     labels=levels(kk$time)) +
  theme_bw() +
  theme(legend.position="bottom")

f <- file.path(CD4CHIC.OUT,"paper/extended-timecourse.png")
ggsave(f,height=8,width=7)
system(paste("display ",f))


head(RESULTS)
head(short)


library(data.table)
library(magrittr)
m <- as.data.table(x)
m[,probe:=rownames(x)]
m <- melt(m,"probe")


m <- merge(m,p,by="variable")
m0 <- m[time=="0h",]
m1 <- m[time!="0h",]
m <- merge(m1,m0[,.(sample.id,probe,value)],by=c("sample.id","probe"))
m[,logFC:=value.x - value.y, by=c("probe","cell","time")]
#m <- unique(m,by=c("probe","cell","time"))

## look up genes in our modules
source("common.R")
modules <- get.modules()

xtract <- function(mod) {
    g <- modules[module==mod, ]$geneSymbol %>% as.character()
    probes <- f[GENE_SYMBOL %in% g, ]$ID
    tmp <- m[probe %in% probes,]
    tmp[,module:=mod]
    tmp <- tmp[,avg.logFC:=mean(logFC), by=c("cell","time","sample.id")]
    ## unique(tmp,by=c("cell","time"))
    tmp
}

mods <- unique(modules$module)
## tscale <- function(x) t(x) %>% scale() %>% t()
## X1 <- lapply(mods, function(mod) tscale(x1) %>% xtract(.,f1,mod))
## X2 <- lapply(mods, function(mod) tscale(x2) %>% xtract(.,f2,mod))
    f <- as.data.table(f)
M <- lapply(mods, function(mod) xtract(mod)) %>% do.call("rbind",.)
M[,mind:=factor(paste(sample.id,module))]
M0 <- unique(M,by=c("cell","module","sample.id"))
M0$avg.logFC <- 0
M0$time <- "0h"
M <- rbind(M,M0)
M <- M[order(M$time),]
Mu <- unique(M,by=c("time","module","cell","sample.id"))
Mu <- Mu[order(Mu$time),]

library(ggplot2)
library(cowplot)
colScale <- make.colScale(unique(M$module))
ggplot(Mu,
       aes(x=time,y=avg.logFC,col=module,group=mind)) + geom_hline(yintercept=0,linetype="dashed") + geom_point() + geom_path() + facet_wrap(~cell ) + colScale + theme_bw()

ggsave(file.path(CD4CHIC.OUT,"paper/extended-timecourse-by-module.pdf"),height=8,width=8)
ggsave(file.path(CD4CHIC.OUT,"paper/extended-timecourse-by-module.tiff"),height=8,width=8)
write.table(Mu,file=file.path(CD4CHIC.OUT,"paper/extended-timecourse-by-module.csv"))
## probes <- f[GENE_SYMBOL=="IL2RA",]$ID
## M[,sind:=paste(sample.id,probe)]
## ggplot(M[probe %in% probes,],
##        aes(x=time,y=value.y - value.x,col=module,group=sind)) + geom_point() + geom_path() + colScale + facet_wrap(sample.id ~ cell)


## library(pheatmap)
## Mo <- unique(M,by=c("time","cell","probe"))
## kk <- unique(Mo[time!="0h",], by=c("cell","time","probe"))
##     i <- as.numeric(as.factor(kk$probe))
##     j <- as.numeric(as.factor(paste(kk$cell,kk$time)))
## mat <- matrix(NA,max(i),max(j))
##    mat[cbind(i,j)] <- kk$logFC
## dimnames(mat) <- list(levels(as.factor(kk$probe)),
##                       levels(as.factor(paste(kk$cell,kk$time))))
## rows <- as.data.frame(unique(Mo[,.(probe,module)]))
## rownames(rows) <- rows$probe
## rows <- rows[rownames(mat),"module",drop=FALSE]
## cols <- as.data.frame(unique(Mo[time!="0h",.(cell,time)]))
## rownames(cols) <- paste(cols$cell,cols$time)
## cols <- cols[colnames(mat),]

## use <- which(rows$module!="grey")
## use <- split(use,rows$module[use]) %>% unlist()
## cuse <- cols$cell %in% c("Th1","Th2")

## pheatmap(mat[use,cuse],annotation_row=rows[use,"module",drop=FALSE],
##          annotation_col=cols,show_rownames=FALSE)



