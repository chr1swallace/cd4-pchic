library(DESeq2)
library(ggplot2)
library(reshape2)
library(plyr)

files <- c("data_gne.202.csv", "data_gne.201.csv", "data_gne.204.csv", "data_gne.203.csv")
condition <- factor(c("nonactivated", "activated", "nonactivated", "activated"), levels=c("nonactivated", "activated"))
pool <- factor(c("1","1","2","2"))
name <- mdply(data.frame(a=condition, b=pool), 'paste', sep='-')[,3]
sampleTable <- data.frame(sampleName=1:4, fileName=files, condition=condition, pool=pool, name=name, stringsAsFactors=F)

fct <- function(file, data, ytext, ltext) {
#  mal <- cbind(data, sampleTable[c(2,1,4,3),])
st <- data.frame(name=factor(c("activated-1", "nonactivated-1", "activated-2", "nonactivated-2"), levels=c("nonactivated-1","activated-1","nonactivated-2","activated-2")))
mal <- cbind(data, st)
  mml <- melt(mal, id.vars=colnames(st))
  pdf(file)
  print(ggplot(mml, aes(x=name, y=value, fill=variable)) + xlab('Sample') + ylab(ytext) + guides(fill=guide_legend(title=ltext)) +
		  geom_bar(stat="identity", position="fill") + theme(axis.text.x = element_text(angle = 90)))
  dev.off()
}

aln <- read.csv('data_aln.csv')
fct('data_anl.map.pdf', as.data.frame(data.matrix(aln)[,1] / 100 * data.matrix(aln)[,2:7]), 'Read ratio', 'Mapping class')
fct('data_anl.ftc.pdf', read.csv('data_ftc.csv'), 'Read count', 'Feature class')

ftc <- c("data_ftc.202.csv", "data_ftc.201.csv", "data_ftc.204.csv", "data_ftc.203.csv")
cnt <- do.call(cbind,lapply(ftc, function(x) read.table(x, header=T)))[, c(1:4 * 7)]
colnames(cnt) <- sampleTable$name
png('data_anl.crr.png', width=2048, height=2048)
pairs(log2(as.matrix(cnt)+0.5), pch=20)
dev.off()

dds <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=".", design= ~ pool + condition)
dds <- estimateSizeFactors(dds)
rld <- rlogTransformation(dds, blind=T)

pdf('data_gdf.pca.pdf')
pca <- plotPCA(rld, intgroup=c("condition", "pool"), returnData=T)
percentVar <- round(100*attr(pca,"percentVar"))
ggplot(pca, aes(PC1, PC2, color=condition, shape=pool)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

dds <- estimateDispersions(dds)

pdf('data_gdf.dsp.pdf')
plotDispEsts(dds)
dev.off()

dds <- nbinomWaldTest(dds)
res <- results(dds)

ngc <- as.data.frame(counts(dds, normalized=T))
colnames(ngc) <- sampleTable$name
ens <- cbind(read.table('data_ref.ens', header=T), as.data.frame(res), ngc)
write.table(ens[order(ens$padj),], file='data_gdf.csv', quote=F, sep='\t', row.names=F)

pdf('data_gdf.maa.pdf')
plotMA(dds,ylim=c(-2,2),main="DESeq2")
dev.off()

pdf('data_gdf.htm.pdf')
top <- ens[order(ens$padj), c("gene", sampleTable$name)][1:20,]
top$gene <- factor(as.vector(top$gene), levels=ens[order(ens$padj),"gene"][1:20])
ggplot(melt(top, id.vars="gene"), aes(y=gene, x=variable)) + geom_tile(aes(fill=value)) + 
		 scale_fill_gradient(name="read count", low="green", high="red", trans="log", labels=function(x) as.character(signif(x,2))) + 
		 ylim(rev(levels(top$gene))) + xlab("sample")
dev.off()
