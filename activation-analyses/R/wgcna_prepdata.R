library(WGCNA)
library(org.Hs.eg.db)
library(data.table)

x <- org.Hs.egENSEMBL
xx <- as.list(org.Hs.egENSEMBL2EG)

options(stringsAsFactors = FALSE);
enableWGCNAThreads(20)
setwd(file.path(CD4CHIC.ROOT,"activation-analyses/R"))
source("common.R")

## get rna-seq limma data @ 4hrs
base.dir<-'/home/oliver/STATS/PROJECTS/T-CELL_PROFILE/WGCNA/'
support.dir<-paste0(base.dir,'support/')

## limma.thresh = 0.05 # SHOULD THIS BE 0.01?
## cd4.limma<-fread(paste0(support.dir,"cd4-limma.csv"),header=TRUE,stringsAsFactors=FALSE,sep=" ")
##
## list of protein_coding genes that are upregulated
##cd4.limma<-subset(cd4.limma,biotype=="protein_coding")
## cd4.down<-subset(cd4.limma, adj.P.Val<=limma.thresh & logFC>0)$id
## cd4.up<-subset(cd4.limma, adj.P.Val<=limma.thresh & logFC<0)$id
## cd4.not.detected<-subset(cd4.limma, is.na(adj.P.Val))$id
## cd4.other<-setdiff(cd4.limma$id, c(cd4.down,cd4.up,cd4.not.detected))

## setattr(cd4.limma, 'colnames', make.names(colnames(cd4.limma)))

## library(ggplot2)
## cd4.limma <- within(cd4.limma,
##                     {nonact.total <- nonact.1+nonact.2
##                      act.total <- act.1+act.2})
## library(ggplot2)
## with(cd4.limma, qplot(-logFC,-log10(adj.P.Val),col=log(Length))) + labs(x="Fold change (act/nonact)")
## with(cd4.limma, qplot(AveExpr,-log10(adj.P.Val))) + geom_smooth()
## with(cd4.limma, qplot(log(nonact.total+1),-log10(adj.P.Val))) + geom_smooth()
## with(cd4.limma, qplot(log(act.total+1),-log10(adj.P.Val))) + geom_smooth()
## with(cd4.limma, qplot(log(nonact.total+1),log(act.total+1),col=-log10(adj.P.Val)))

## load the expression data and meta data
tcell.expr.data.file<-paste0(support.dir,'/t-cell-act.profiles.RData')
(load(file=tcell.expr.data.file))
## load the meta data
tcell.expr.meta.data.file<-paste0(support.dir,'/Xaq_PTPN22_pairs_genotype_status-1.csv')
pheno<-read.csv(file=tcell.expr.meta.data.file,header=TRUE,stringsAsFactors=FALSE)
pheno$index<-with(pheno,paste(Pair,Donor,sep="."))


## get the log2 transformed data
exp.dataset<-exprs(gene.vsnorm.data)

## get a list of samples as we want to transform
samples<-pData(gene.vsnorm.data)


## pull these apart to add meta data information more complicated
## as id's don't follow a strict pattern
sample.df<-do.call("rbind",lapply(strsplit(rownames(samples),"\\_"),function(x){
	sname<-paste(x,sep="_",collapse="_")
	x[1]<-sub("Pair","",x[1])
	x[2]<-sub("t","",x[2])
	if(length(x)==4)
		return(data.frame(sname=sname,pair=x[1],time=x[2],treatment=x[3],individual=sub("^D([^\\.]+.*)\\.CEL","\\1",x[4])))
	if(length(x)==3)
		return(data.frame(sname=sname,pair=x[1],time=x[2],treatment='T0',individual=sub("^D([^\\.]+.*)\\.CEL","\\1",x[3])))
	#x[3]<-paste(x[3],x[4],sep="_",collapse="_")
	return(data.frame(sname=sname,pair=x[1],time=x[2],treatment=x[3],individual=sub("^D([^\\.]+.*)\\.CEL","\\1",x[4])))
}))

##fixes D prefix inclusion
sample.df$individual<-sub("D","",sample.df$individual)
##fixed instances where t0 is not named as US
sample.df[sample.df$time=="0" & sample.df$treatment=="T0",]$treatment="US"
sample.df$index<-with(sample.df,paste(pair,individual,sep="."))
pheno.data<-merge(pheno,sample.df,by.x="index",by.y="index",all.x=TRUE)

## expression set filtering
genes.to.remove<-which(apply(exp.dataset,1,max)<expression.threshold)
if(length(genes.to.remove)>0){
	filt.exp<-exp.dataset[-genes.to.remove,]
}else{
	filt.exp<-exp.dataset
}

failed.hybs <-paste0(c("Pair15_t0_D2","Pair15_t2_S_D1","Pair15_t2_S_D2","Pair15_t4_S_D1"),".CEL")
failed.hybs %in% colnames(exp.dataset)
filt.exp <- filt.exp[, setdiff(colnames(exp.dataset),failed.hybs)]

## for the time being remove any probes that don't have an ensgene id
## array.annotation.file
array.annot.file<-paste0(support.dir,'HuGene-1_0-st-v1.e67.genes2GN-prbsts.tab')
map=read.table(array.annot.file,header=TRUE,stringsAsFactors=FALSE)
ens67.biotype.file<-paste0(support.dir,'e67_gene_biotyes.csv')
biot<-read.csv(ens67.biotype.file)
map<-merge(map,biot,by.x='ens.gene.id',by.y='Ensembl.Gene.ID',all.x=TRUE)
map.filt<-subset(map,!is.na(ens.gene.id) & Gene.Biotype=="protein_coding")
## remove those ensid's that map to multiple entrez id's

## for GO enrichment WGCNA need entrez id code below 
## makes sure that this all works at the cost of removing
## genes
#map.filt$entrezId<-xx[map.filt$ens.gene.id]
#map.filt<-map.filt[sapply(map.filt$entrezId,length)==1,]
#map.filt$entrezId<-unlist(map.filt$entrezId)
## removed duplicated entrez id's (prob need to do this better in future)
#map.filt<-map.filt[!duplicated(map.filt$entrezId),]
filt.exp<-filt.exp[which(rownames(filt.exp) %in% map.filt$probeset.id),]


## lookup for relating pheno to colnames -- here position in matrix is location in pheno.data
## value is column in filt.exp

## expset2pheno.data<-match(pheno.data$sname,colnames(filt.exp))
## pheno.data$expsetmatch<-expset2pheno.data
pheno.data <- pheno.data[ match(colnames(filt.exp), pheno.data$sname), ]

with(pheno.data, ftable(time,treatment))
with(pheno.data, ftable(time,treatment)) %>% sum()

## compute change from baseline for each sample
## used for plotting resultant modules
## DE matrix
samples.t0 <- subset(pheno.data,time==0)
samples.t1 <- subset(pheno.data,time!=0 & uniqueID %in% samples.t0$uniqueID)
samples.t1 <- merge(samples.t1,samples.t0[,c("uniqueID","sname")],by="uniqueID",suffixes=c("",".0"))

all(samples.t1$sname %in% colnames(filt.exp))
all(samples.t1$sname.0 %in% colnames(filt.exp))
treated.t1 <- subset(samples.t1,treatment=="S")

delta.exp <- filt.exp[,samples.t1$sname] - filt.exp[,samples.t1$sname.0] 
delta.treated <- filt.exp[,treated.t1$sname] - filt.exp[,treated.t1$sname.0] 

## switch this in to work with changes in expression wrt baseline (t=0 unstimulated)
s.exp<-list(all=t(delta.exp),
            treated=t(delta.treated))

## correct sample mix up


pca <- prcomp(t(filt.exp))
tmp <- pheno.data
tmp$pc1 <- pca$x[,"PC1"]
tmp$pc2 <- pca$x[,"PC2"]
ggplot(tmp,aes(x=pc1,y=pc2,col=treatment,pch=time)) + geom_point()
subset(tmp, pc2 > -10 & time==21 & treatment=="S")
subset(tmp, pc2 < -20 & time==21 & treatment=="US")

pheno.data[pheno.data$sname %in% c("Pair10_t21_S_D1.CEL","Pair10_t21_S_D2.CEL"),"treatment"] <- "US"
pheno.data[pheno.data$sname %in% c("Pair10_t21_US_D1.CEL","Pair10_t21_US_D2.CEL"),"treatment"] <- "S"

save(filt.exp,s.exp,map.filt,pheno.data,file=file.expression)
