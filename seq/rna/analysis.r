library(ggplot2)
library(reshape2)
library(plyr)

design <- data.frame(n=c(202, 201, 204, 203), sample=c('1', '1', '2', '2'), condition=c('nonactivated', 'activated', 'nonactivated', 'activated'), stringsAsFactors=T)
name <- mdply(data.frame(a=design$sample, b=design$condition), 'paste', sep='-')[,3]
design <- cbind(design, data.frame(name=name))

fct <- function(file, data, type, ytext, ltext) {
  mal <- cbind(data, design)
  mml <- melt(mal, id.vars=colnames(design))
  pdf(file)
  print(ggplot(mml, aes(x=name, y=value, fill=variable)) + xlab('Sample') + ylab(ytext) + guides(fill=guide_legend(title=ltext)) +
		  geom_bar(stat="identity", position=type))
  dev.off()
}

alr <- read.csv('data_aln.csv')
aln <- as.data.frame(data.matrix(alr)[,1] * data.matrix(alr)[,2:7])
fct('data_anl.map.rel.pdf', aln, 'fill', 'Read ratio', 'Mapping class')
fct('data_anl.map.abs.pdf', aln, 'stack', 'Read count', 'Mapping class')

ftc <- read.csv('data_ftc.csv')
fct('data_anl.ftc.rel.pdf', ftc, 'fill', 'Read ratio', 'Feature class')
fct('data_anl.ftc.abs.pdf', ftc, 'stack', 'Read count', 'Feature class')

library(Rsubread)
library(limma)
library(edgeR)

alf <- paste('data_aln.', design$n, '.Aligned.out.bam', sep='')
fc <- featureCounts(files=alf, nthreads=16, isPairedEnd=T, requireBothEndsMapped=T, strandSpecific=2, minMQS=40, countChimericFragments=F, 
      annot.ext='/dunwich/scratch/ar756/.local/ensembl/75/Homo_sapiens.GRCh37.75.gtf', isGTFAnnotationFile=T)
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])
isexpr <- rowSums(cpm(x) > 10) >= 2
x <- x[isexpr,]
d <- model.matrix(~design$condition + design$sample)
y <- voom(x,d,plot=TRUE)
pdf('data_anl.mds.pdf')
plotMDS(y)
dev.off()
fit <- eBayes(lmFit(y,d))
topTable(fit,coef=2)
fc$counts[topTable(fit, coef=2, number=25)$GeneID,]

c<- merge(ens, fc$counts, by.x='id', by.y='row.names')
t<-merge(c, topTable(fit, coef=2, number=nrow(ens)), all=T, by.x='id', by.y='GeneID')
