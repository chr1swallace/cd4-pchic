## DONE
## vasculitis, SLE https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-145
## Graves https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-71957
## RA https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-20098 https://www.ncbi.nlm.nih.gov/pubmed/22532634
## Multiple diseases (RNAseq) https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-60424

## TODO
## SLE https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2713/?query=cd4%20AND%20patients&sortby=assays&sortorder=descending

## http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0140049 - cd - 15 vs 11, stim and non-stim - rnaseq

## http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13732 Longitudinal, 37 CIS (MS) vs 28 Controls - but naive CD4
## http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55447 Multiethnic SLE 21 african american,21 European american 5 AA controls and 5 EA controls.

## http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60236 - Interesting and potentially relevant to Burren et al. - I have a feeling we have discussed before - eqtl, not disease



## Jimmy's (UC, CRO)
(load("~/scratch/cd4chic/jimmy-microarray.RData"))
IBD.de <- results[[1]]

library(ArrayExpress)
a1 <- getAE("E-MEXP-145",type="processed",path="~/scratch/cd4chic/disease-expr")
a2 <- getAE("GSE71956",type="processed",path="~/scratch/cd4chic/disease-expr")
a3 <- getAE("E-GEOD-20098",type="processed",path="~/scratch/cd4chic/disease-expr")
a4 <- getAE("E-GEOD-60424",type="processed",path="~/scratch/cd4chic/disease-expr")

reader <- function(a) {
    files <- grep("K3K",a$processedFiles,value=TRUE,invert=TRUE)
    DATA <- vector("list",length(files))
    names(DATA) <- files
    for(i in seq_along(files)) {        
        f <- file.path(a$path, files[i])
        message(f)
        d <- read.table(f,skip=2,header=FALSE)
        h <- scan(f,what="",nlines=1,sep="\t")
        colnames(d) <- h
        DATA[[i]] <- d
    }
    return(DATA)
}

d1 <- reader(a1)
d2 <- reader(a2)
d3 <- reader(a3)
dim(d1[[1]]) # 25 - 8 subjects - almost surely too small


## data 2 - graves
library(magrittr)
length(d2)
ids <- lapply(d2, "[[",1)
ok <- TRUE
for(i in 2:length(ids)) {
    if(!(identical(ids[[i]],ids[[1]])))
        stop("mismatch at ",i)
}
expr <- lapply(d2, "[[", 2) %>% do.call("cbind",.)
rownames(expr) <- ids[[1]]

s2 <- read.table(file.path(a2$path,a2$sdrf),sep="\t",comment.char="",quote="",header=TRUE)
rownames(s2) <- s2$Derived.Array.Data.File
s2 <- s2[colnames(expr),]
use <- s2$FactorValue..cell.type.=="CD4 T cells"

library(ggplot2)
p <- prcomp(t(expr))
qplot(p$x[,1],p$x[,2],col=s2$FactorValue..cell.type.,pch=s2$Characteristics..diagonsis.,geom="point")

e2 <- vsn2(expr)


library(limma)
s <- s2[use,]
design <- with(s, model.matrix(~0 + Characteristics..diagonsis.))
head(design)
colnames(design) <- c("graves","hv")
fit <- lmFit(e2[,use], design)
contrast.matrix <- makeContrasts(graves-hv,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
graves.diff <- topTable(fit2, number=Inf)

f <- file.path(a2$path,a2$adf)
system(paste("head -n20 ",f," | nl"))
probes <- read.table(f,sep="\t",comment.char="",quote="",header=TRUE,skip=14,row.names=1)
probes <- probes[rownames(graves.diff),]
head(probes)
graves.diff <- cbind(graves.diff,probes)

## data 3 - multiple 
length(d3)
ids <- lapply(d3, "[[",1)
ok <- TRUE
for(i in 2:length(ids)) {
    if(!(identical(ids[[i]],ids[[1]])))
        stop("mismatch at ",i)
}
expr <- lapply(d3, "[[", 2) %>% do.call("cbind",.)
rownames(expr) <- ids[[1]]

s3 <- read.table(file.path(a3$path,a3$sdrf),sep="\t",comment.char="",quote="",header=TRUE)
rownames(s3) <- s3$Derived.Array.Data.File

s3 <- s3[colnames(expr),]
head(s3)
with(s3,table(Characteristics.cell.type.)) ## all CD4
with(s3,table(Comment..Sample_source_name.))

library(ggplot2)
p <- prcomp(t(expr))
qplot(p$x[,1],p$x[,2],col=s3$Comment..Sample_source_name.,geom="point")

                                        #e3 <- vsn2(t(expr))
e3 <- expr

s3$pc1 <- p$x[,1]

s3$test1 <- ifelse(s3$Comment..Sample_source_name.=="ex-vivo peripheral blood CD4 T-cells, non-RA/non-inflammatory arthritis","ctl",
            ifelse(grepl("ACPA",s3$Comment..Sample_source_name.),"RA",NA))
s3$test2 <- ifelse(grepl("ACPA",s3$Comment..Sample_source_name.),"RA","ctl")

library(limma)
s <- s3
design1 <- with(s, model.matrix(~0 + pc1 + test1))
design2 <- with(s, model.matrix(~0 + pc1 + test2))
head(design1)
head(design2)
colnames(design1) <- colnames(design2) <- c("pc1","ctl","RA")
fit1 <- lmFit(e3[,!is.na(s3$test1)], design1)
fit2 <- lmFit(e3, design2)
contrast1 <- makeContrasts(RA-ctl,levels=design1)
contrast2 <- makeContrasts(RA-ctl,levels=design2)
fit1 <- eBayes(contrasts.fit(fit1, contrast1))
fit2 <- eBayes(contrasts.fit(fit2, contrast1))
RA.diff1 <- topTable(fit1, number=Inf)
RA.diff2 <- topTable(fit2, number=Inf)

f <- file.path(a3$path,a3$adf)
system(paste("head -n20 ",f," | nl"))
probes <- read.table(f,sep="\t",comment.char="",quote="",header=TRUE,skip=18,row.names=1)
probes <- probes[rownames(RA.diff1),]
head(probes)
RA.diff1 <- cbind(RA.diff1,probes)

probes <- probes[rownames(RA.diff2),]
head(probes)
RA.diff2 <- cbind(RA.diff2,probes)

head(RA.diff1)
head(RA.diff2)

library(ggplot2)
plotter <- function(x) ggplot(x[x$adj.P.Val<0.1,],aes(x=logFC,y=-log10(adj.P.Val))) + geom_point()

plotter(RA.diff1)
plotter(RA.diff2)
plotter(graves.diff)

RA.diff <- RA.diff2
save(RA.diff,file="~/gscratch/RA-CD4-diff.RData")

 ##    IGFL2
 ##      GNG4
 ##     SOCS3
 ##    CXCL13
 ## LOH11CR2A
 ##       EDA
 ##      PIM1
 ##     RNU64
 ##     SOCS3
 ##     CCNE2
 ##    ACVR2A
 ##   SNORA10
 ##     PDCD1
